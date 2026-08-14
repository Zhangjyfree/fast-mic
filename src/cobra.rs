//! FBA + CycleFreeFlux + pFBA, with lexicographic max-min co-culture.
//!
//! Design goals (vs. previous version):
//!   1. NO empirical biological hyperparameters in the algorithm.
//!      All remaining numerical constants are LP solver tolerances.
//!   2. Energy-cycle removal is mathematically principled (CFF, Desouki
//!      et al. 2015) and applied unconditionally: it is a no-op on
//!      cycle-free solutions, so it cannot hurt correct answers.
//!   3. Co-culture allocation uses lexicographic max-min (Rawlsian
//!      egalitarian, then utilitarian), which has a unique optimum
//!      and no scan / tradeoff search.
//!   4. Cross-feeding analysis uses sign-only rules instead of
//!      magnitude thresholds.
//!   5. Each emitted MetaboliteExchange carries an `inferred` flag
//!      separating direct (mass-balance, sign-based) from facilitation
//!      (benefit-driven fallback) attributions, so downstream tools
//!      can filter on evidence quality without re-deriving the rule.
//!
//! ## Locked-CFF entry point (`run_fba_locked`)
//!
//! `run_fba` is the legacy entry. `run_fba_locked` extends it with an
//! `extra_locks: &[String]` parameter that *additionally* fixes a list
//! of reaction fluxes at their FBA-pass values during the CFF / pFBA
//! pass.  In the pairwise co-culture context this is used to pin each
//! species' biomass to its phase-1 value while CFF removes interior
//! cycles in the merged model.

use crate::medium;
use crate::model::{
    self, InteractionType, MetabolicModel, Metabolite, MetaboliteExchange,
    PairwiseResult, Reaction, SharedResource, TargetFluxRecord,
};
use good_lp::{variable, ProblemVariables, Expression, SolverModel, Solution, Variable, highs};
use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicUsize, Ordering};

// ============================================================
// Configurable parameters
// ============================================================

#[derive(Debug, Clone)]
pub struct AnalysisParams {
    /// Tolerance for pinning reaction fluxes to FBA-pass values in CFF /
    /// pFBA lock constraints: `v ∈ [v* − lock_tol, v* + lock_tol]`.
    ///
    /// Default 1e-5 (10× `NUMERICAL_TOL`, 100× the HiGHS feasibility
    /// tolerance 1e-7).  Drop to 1e-7 only when exact reproduction of
    /// single-species FBA values is required.
    pub lock_tol: f64,
    /// Co-culture growth-allocation objective.  `false` (default) uses the
    /// lexicographic max-min (Rawlsian + utilitarian) allocation; `true`
    /// pins the co-culture growth ratio to the monoculture ratio
    /// (μ_A/μ_B = μ_A^alone/μ_B^alone) for the objective-sensitivity analysis.
    pub fixed_ratio: bool,
}

impl Default for AnalysisParams {
    fn default() -> Self {
        Self {
            lock_tol: LOCK_TOL,
            fixed_ratio: false,
        }
    }
}

// ============================================================
// Numerical tolerances (LP solver-driven, NOT biology)
// ============================================================

/// Biological viability floor (h⁻¹). ≈ 4 doublings/day; well above the
/// solver feasibility tolerance.
const MIN_VIABLE_GROWTH: f64 = 1e-4;

/// Single numerical tolerance used everywhere LP outputs are compared:
/// flux-zero detection, bound-saturation tests, and fairness-floor checks.
/// Set 10× the HiGHS default primal feasibility tolerance (1×10⁻⁷;
/// Huangfu & Hall 2018) to provide a safety margin against solver round-off.
const NUMERICAL_TOL: f64 = 1e-6;

/// Tolerance used when pinning reaction fluxes to their FBA-pass values in
/// CycleFreeFlux / pFBA lock constraints:
///   v ∈ [v* − LOCK_TOL, v* + LOCK_TOL]
///
/// Set 10× NUMERICAL_TOL (= 100× the HiGHS feasibility tolerance).  Drop
/// to 1e-7 only when bit-exact reproduction of single-species FBA values
/// is required.
const LOCK_TOL: f64 = 1e-5;

// ============================================================
// FBA Result
// ============================================================

#[derive(Debug, Clone)]
pub struct FBAResult {
    pub objective_value: f64,
    pub flux_distribution: HashMap<String, f64>,
    /// Stage-1 FBA-optimal biomass, carried through the CycleFreeFlux/pFBA
    /// step. Equals `objective_value` when no CFF refinement was applied;
    /// otherwise `fba_optimal - objective_value` is the post-CFF biomass
    /// deviation reported as a loop-removal validation metric.
    /// Used by the `bench-cff-deviation` binary, so silence dead-code warnings
    /// when compiling the main `fast-mic` binary.
    #[allow(dead_code)]
    pub fba_optimal: f64,
}

// ============================================================
// Sparse stoichiometric matrix
// ============================================================

#[derive(Debug, Clone)]
pub struct SparseS {
    pub rows: Vec<Vec<(usize, f64)>>,
}

impl SparseS {
    pub fn build(reactions: &[Reaction], metabolites: &[Metabolite]) -> Self {
        let met_index: HashMap<&str, usize> = metabolites
            .iter().enumerate().map(|(i, m)| (m.id.as_str(), i)).collect();
        let mut rows: Vec<Vec<(usize, f64)>> =
            (0..metabolites.len()).map(|_| Vec::new()).collect();
        for (j, rxn) in reactions.iter().enumerate() {
            for (met_id, coeff) in &rxn.metabolites {
                if let Some(&i) = met_index.get(met_id.as_str()) {
                    rows[i].push((j, *coeff));
                }
            }
        }
        for row in &mut rows {
            if row.len() < 2 { continue; }
            row.sort_by_key(|(j, _)| *j);
            let mut merged: Vec<(usize, f64)> = Vec::with_capacity(row.len());
            for &(j, c) in row.iter() {
                if let Some(last) = merged.last_mut() {
                    if last.0 == j { last.1 += c; continue; }
                }
                merged.push((j, c));
            }
            *row = merged;
        }
        SparseS { rows }
    }
}

fn add_mass_balance_constraints(
    mut lp: good_lp::solvers::highs::HighsProblem,
    sparse_s: &SparseS, flux_vars: &[Variable],
) -> good_lp::solvers::highs::HighsProblem {
    for row in &sparse_s.rows {
        if row.is_empty() { continue; }
        let expr: Expression = row.iter().map(|&(j, c)| c * flux_vars[j]).sum();
        lp = lp.with(expr.eq(0.0));
    }
    lp
}

// ============================================================
// Generic helpers
// ============================================================

#[inline]
fn snap_zero(x: f64) -> f64 { if x.abs() < 1e-9 { 0.0 } else { x } }

/// Prefix-only exchange reaction check for use in CFF/pFBA locking and
/// cross-feeding attribution, where reaction names follow the `EX_*` /
/// `R_EX_*` convention by construction (merged pairwise co-culture models).
/// For canonical COBRApy-style checking use `model::is_canonical_exchange`.
#[inline]
pub fn is_exchange_id(id: &str) -> bool {
    id.starts_with("R_EX_") || id.starts_with("EX_")
}

#[inline]
fn strip_compartment_suffix(s: &str) -> &str {
    // Covers both ModelSEED (`_e0`, `_c0`, `_p0`) and BiGG-style
    // (`_e`, `_c`, `_p`, `_m`) compartments, plus legacy COBRA
    // bracket form (`[e]`, `[c]`, `[p]`, `[m]`).
    for suf in &["_e0", "_c0", "_p0", "_e", "_c", "_p", "_m",
                 "[e]", "[c]", "[p]", "[m]"] {
        if let Some(stripped) = s.strip_suffix(suf) { return stripped; }
    }
    s
}

/// Canonical metabolite base ID: strip `M_` prefix, lowercase, drop compartment suffix.
/// Use this to look up metabolites in any of the filtering HashSets defined below
/// (`inorganic_set`, `fba_artifact_set`).
fn metabolite_base(id: &str) -> String {
    let clean = id.strip_prefix("M_").unwrap_or(id).to_lowercase();
    strip_compartment_suffix(&clean).to_string()
}


fn has_real_gene_support(genes: &[String]) -> bool {
    if genes.is_empty() { return false; }
    genes.iter().any(|g| {
        let lower = g.to_lowercase();
        !lower.contains("spontaneous") && !lower.starts_with("s0001") && lower != "unknown"
    })
}

fn compute_gene_supported_fraction(
    cf_score: f64,
    a_to_b: &[MetaboliteExchange],
    b_to_a: &[MetaboliteExchange],
) -> f64 {
    if cf_score < NUMERICAL_TOL { return f64::NAN; }
    let all: Vec<&MetaboliteExchange> = a_to_b.iter().chain(b_to_a.iter()).collect();
    if all.is_empty() { return f64::NAN; }
    let total_flux: f64 = all.iter().map(|e| e.flux).sum();
    if total_flux <= NUMERICAL_TOL { return f64::NAN; }
    let supported_flux: f64 = all.iter()
        .filter(|e| has_real_gene_support(&e.donor_genes)
            || has_real_gene_support(&e.receiver_genes))
        .map(|e| e.flux).sum();
    (supported_flux / total_flux).clamp(0.0, 1.0)
}

// ============================================================
// LP variable creation
// ============================================================

fn create_flux_var(vars: &mut ProblemVariables, r: &Reaction, _params: &AnalysisParams) -> Variable {
    vars.add(variable().min(r.lower_bound).max(r.upper_bound))
}

fn create_coculture_flux_var(
    vars: &mut ProblemVariables, r: &Reaction,
    a_viable: bool, b_viable: bool, _params: &AnalysisParams,
) -> Variable {
    let is_dead_a = !a_viable && r.id.ends_with("_A");
    let is_dead_b = !b_viable && r.id.ends_with("_B");
    if is_dead_a || is_dead_b { return vars.add(variable().min(0.0).max(0.0)); }
    vars.add(variable().min(r.lower_bound).max(r.upper_bound))
}

fn sorted_reactions(model: &MetabolicModel) -> Vec<Reaction> {
    let mut r: Vec<Reaction> = model.reactions.values().cloned().collect();
    r.sort_by(|a, b| a.id.cmp(&b.id)); r
}

fn sorted_metabolites(model: &MetabolicModel) -> Vec<Metabolite> {
    let mut m: Vec<Metabolite> = model.metabolites.values().cloned().collect();
    m.sort_by(|a, b| a.id.cmp(&b.id)); m
}

fn inorganic_set() -> HashSet<String> {
    [
        // ── BiGG-style names ──────────────────────────────────────────────
        "h2o","h","h2","o2","co2",
        "nh4","pi","so4",
        "h2s","hco3",
        "mg2","ca2","k","na1",
        "fe2","fe3","zn2","mn2","cu2","cobalt2",
        "cl","mobd","ni2","sel","tungs","slnt",
        // ── ModelSEED cpd IDs (these were the gap that let H2O/Pi/H+/NH4
        //    leak into >17 000 cross-feeding rows in the gut output) ───────
        "cpd00001",  // H2O
        "cpd00067",  // H+ (proton)
        "cpd00007",  // O2
        "cpd00011",  // CO2
        "cpd00013",  // NH4+ (ammonium)
        "cpd00009",  // Pi (orthophosphate)
        "cpd00048",  // SO4 (sulfate)
        "cpd00239",  // H2S (hydrogen sulfide)
        "cpd00242",  // HCO3- (bicarbonate)
        "cpd00099",  // Cl-
        "cpd00254",  // Mg2+
        "cpd00063",  // Ca2+
        "cpd00205",  // K+
        "cpd00971",  // Na+
        "cpd10515",  // Fe2+
        "cpd10516",  // Fe3+
        "cpd00034",  // Zn2+
        "cpd00030",  // Mn2+
        "cpd00058",  // Cu2+
        "cpd00149",  // Co2+
        "cpd00244",  // Ni2+
        "cpd11574",  // Molybdate (MoO4)
        "cpd00528",  // Nitrogen gas (N2)
        "cpd00209",  // Nitrate (NO3-)
        "cpd00075",  // Nitrite (NO2-)
        "cpd00531",  // Hg2+
        "cpd00264",  // Selenate
        "cpd00048",  // already listed (sulfate)
        // Added after empirical review of v1 results: PPi dominated several
        // Akk-import panels but is a textbook FBA artifact (cytoplasmic
        // pyrophosphatase hydrolyses it to 2 Pi instantly; no transporter).
        "cpd00012",  // PPi (pyrophosphate)
        "ppi",       // BiGG alias
    ].iter().map(|s| s.to_string()).collect()
}

// ── Filter counters (process-global, accumulated across all pair calls) ─────
// Exposed via `print_filter_summary()` so the binary can report exactly how
// many exchanges were filtered as artifacts vs kept as real cross-feeding.
static FILTER_EVALUATED:       AtomicUsize = AtomicUsize::new(0);
static FILTER_INORGANIC:       AtomicUsize = AtomicUsize::new(0);
static FILTER_COFACTOR:        AtomicUsize = AtomicUsize::new(0);
static FILTER_REACTIVE:        AtomicUsize = AtomicUsize::new(0);
static FILTER_ELECTRON:        AtomicUsize = AtomicUsize::new(0);
static FILTER_COA_DERIV:       AtomicUsize = AtomicUsize::new(0);
static FILTER_BOUND_SATURATED: AtomicUsize = AtomicUsize::new(0);
static FILTER_KEPT:            AtomicUsize = AtomicUsize::new(0);

/// Print a summary of how many candidate exchange metabolites were filtered
/// as artifacts vs kept as real cross-feeding signal.  Call once at the end
/// of the run from `main.rs`.
pub fn print_filter_summary() {
    let evaluated = FILTER_EVALUATED.load(Ordering::Relaxed);
    if evaluated == 0 { return; }
    let inorganic = FILTER_INORGANIC.load(Ordering::Relaxed);
    let cofactor  = FILTER_COFACTOR.load(Ordering::Relaxed);
    let reactive  = FILTER_REACTIVE.load(Ordering::Relaxed);
    let electron  = FILTER_ELECTRON.load(Ordering::Relaxed);
    let coa       = FILTER_COA_DERIV.load(Ordering::Relaxed);
    let bound     = FILTER_BOUND_SATURATED.load(Ordering::Relaxed);
    let kept      = FILTER_KEPT.load(Ordering::Relaxed);
    let artifacts = cofactor + reactive + electron + coa;
    println!("\n🧪 [FILTER SUMMARY]");
    println!("  Total candidate exchanges evaluated: {}", evaluated);
    println!("  Filtered as inorganics (H+, H2O, CO2, NH4, Pi, metals): {}", inorganic);
    println!("  Filtered as FBA artifacts: {}", artifacts);
    println!("    - Cofactors / energy carriers (CoA, ATP, NAD(P), FAD, THF): {}", cofactor);
    println!("    - Acyl-CoA derivatives (name-based, *-CoA): {}", coa);
    println!("    - Reactive intermediates (aldehydes, methylglyoxal, hpyr): {}", reactive);
    println!("    - Electron carriers (menaquinol/-quinone, ferredoxin): {}", electron);
    println!("  Filtered as bound-saturated (medium boundary): {}", bound);
    println!("  Kept as candidate cross-feeding metabolites: {}", kept);
}

/// Metabolites that appear in ModelSEED/CarveMe exchange reactions for
/// mass/charge balancing but are not biologically transferred between cells.
/// Excluding them prevents spurious "cross-feeding" hits in the output.
///
/// Four categories (see methods):
///   1. Cofactors / energy carriers (CoA, ATP/ADP/AMP, NAD(P)(H), FAD, THF series)
///   2. Reactive intermediates / pseudo-products (formaldehyde, acetaldehyde,
///      glyceraldehyde, methylglyoxal, hydroxypyruvate, glycolaldehyde, glyoxylate)
///   3. Technical electron carriers (menaquinol/-quinone, ferredoxin)
///   4. Real cross-feeding metabolites (SCFAs, amino acids, vitamins,
///      nucleosides, sugars, bile acids, indoles) — never added here.
///
/// All entries are stored in *base* form (no `M_` prefix, no `_e/_c/_p/_e0`
/// suffix, lowercase) so a single `metabolite_base()` lookup suffices.
/// Note: acyl-CoA derivatives (~900 in ModelSEED) are caught by the
/// name-based `is_coa_derivative_name()` helper, not enumerated here.
fn cofactor_set() -> HashSet<String> {
    [
        // ModelSEED IDs
        "cpd00010",  // CoA
        "cpd00002",  // ATP
        "cpd00008",  // ADP
        "cpd00018",  // AMP
        "cpd00003",  // NAD
        "cpd00004",  // NADH
        "cpd00006",  // NADP
        "cpd00005",  // NADPH
        "cpd00015",  // FAD
        "cpd00982",  // FADH2
        "cpd00087",  // THF
        "cpd00125",  // 5,10-methylene-THF
        "cpd00201",  // 10-formyl-THF
        "cpd00345",  // 5-methyl-THF
        // BiGG aliases
        "coa", "atp", "adp", "amp", "nad", "nadh", "nadp", "nadph",
        "fad", "fadh2", "thf", "mlthf", "5mthf", "10fthf",
    ].iter().map(|s| s.to_string()).collect()
}

fn reactive_intermediate_set() -> HashSet<String> {
    // NOTE: Glyoxylate (cpd00040 / glx) intentionally NOT filtered -- it's the
    // central metabolite of the glyoxylate cycle and can genuinely accumulate
    // extracellularly in anaerobes / acetogens / Clostridia. Discuss in text.
    [
        "cpd00055",  // Formaldehyde
        "cpd00071",  // Acetaldehyde
        "cpd00448",  // D-Glyceraldehyde
        "cpd00428",  // Methylglyoxal
        "cpd00145",  // Hydroxypyruvate
        "cpd00229",  // Glycolaldehyde
        "acald", "fald", "glyald", "mthgxl", "hpyr", "gcald",
    ].iter().map(|s| s.to_string()).collect()
}

/// "Sometimes-real-sometimes-artifact" cross-feeding metabolites.
/// We DON'T filter these (they can be genuine cross-feeds in some anaerobes)
/// but we flag them as low-confidence so the figures can render them
/// differently and reviewers see we considered the artifact risk.
///
/// Pyruvate: cytosolic intermediate, rarely transported but observed in some
/// Prevotella / Bacteroides co-cultures.
/// Lactate (L/D): real in mixed-acid fermenters; ambiguous elsewhere.
/// Ethanol: real only for Saccharomyces / yeast-bacteria pairings.
/// Glycerol: frequently a "leaky boundary" in anaerobe models.
fn low_confidence_crossfeed_set() -> HashSet<String> {
    [
        "cpd00020", "pyr",        // Pyruvate
        "cpd00159", "lac__l",     // L-Lactate
        "cpd00221", "lac__d",     // D-Lactate
        "cpd00363", "etoh",       // Ethanol
        "cpd00100", "glyc",       // Glycerol
    ].iter().map(|s| s.to_string()).collect()
}

fn electron_carrier_set() -> HashSet<String> {
    [
        "cpd11451",  // Menaquinol
        "cpd11606",  // Menaquinone 7
        "cpd11620",  // Reduced ferredoxin
        "cpd11621",  // Oxidized ferredoxin
        "mqn7", "mqn8", "mql7", "mql8",
        "fdxo", "fdxr", "fdxox", "fdxrd",
    ].iter().map(|s| s.to_string()).collect()
}

/// Union of all three artifact sub-categories.  Convenience for tests and
/// callers who don't need per-category attribution; production code uses
/// the individual sub-sets so the filter counters can attribute correctly.
#[allow(dead_code)]
fn fba_artifact_set() -> HashSet<String> {
    let mut s = cofactor_set();
    s.extend(reactive_intermediate_set());
    s.extend(electron_carrier_set());
    s
}

/// Name-based detection of acyl-CoA derivatives (Acetyl-CoA, Propionyl-CoA,
/// Succinyl-CoA, Malonyl-CoA, ~900 total in ModelSEED).
/// Catches anything whose human-readable name ends in `-CoA` (case-insensitive),
/// while whitelisting metabolites that *contain* "CoA" in their name but are
/// genuine cross-feeding species (CoA precursors, vitamin B5, etc.).
fn is_coa_derivative_name(name: &str) -> bool {
    let lc = name.to_lowercase();

    // Whitelist: CoA biosynthesis precursors and adjacent vitamins that
    // ARE biologically exchanged and must remain in the cross-feeding output.
    let whitelist = [
        "pantothenate",          // Vitamin B5 / pantothenic acid
        "pantothenic",
        "vitamin b5",
        "pantetheine",           // CoA precursor, can be exchanged
        "4'-phosphopantetheine", // CoA biosynthesis intermediate (border case)
        "phosphopantetheine",
    ];
    if whitelist.iter().any(|w| lc.contains(w)) { return false; }

    lc.ends_with("-coa")
        || lc.ends_with(" coa")
        || lc.ends_with("coenzyme a")
        || lc.ends_with("coenzyme-a")
        || lc.starts_with("acyl-coa")
}

/// Find the index of the biomass reaction in `reactions`.
///
/// Returns `None` when the reaction cannot be identified — the caller must
/// bail rather than silently fall back to `reactions[0]`, which is
/// alphabetically the first reaction (usually an exchange or metabolic
/// reaction) and would produce wrong results without any warning.
fn find_biomass_reaction(reactions: &[Reaction], biomass_id: Option<&str>) -> Option<usize> {
    if let Some(bio_id) = biomass_id {
        if let Some(idx) = reactions.iter().position(|r| r.id == bio_id) { return Some(idx); }
    }
    reactions.iter().position(|r| {
        let il = r.id.to_lowercase(); let nl = r.name.to_lowercase();
        il.contains("biomass") || il.contains("growth")
            || nl.contains("biomass") || nl.contains("growth")
    })
}

fn find_merged_biomass_indices(
    reactions: &[Reaction], biomass_a_id: Option<&str>, biomass_b_id: Option<&str>,
) -> anyhow::Result<(usize, usize)> {
    let bio_a = biomass_a_id
        .and_then(|id| reactions.iter().position(|r| r.id == format!("{}_A", id)))
        .or_else(|| reactions.iter().position(|r| {
            let low = r.id.to_lowercase();
            (low.contains("biomass") || low.contains("growth")) && r.id.ends_with("_A")
        }))
        .ok_or_else(|| anyhow::anyhow!(
            "Cannot find biomass reaction for species A (tried id={:?} and name heuristic)",
            biomass_a_id
        ))?;
    let bio_b = biomass_b_id
        .and_then(|id| reactions.iter().position(|r| r.id == format!("{}_B", id)))
        .or_else(|| reactions.iter().position(|r| {
            let low = r.id.to_lowercase();
            (low.contains("biomass") || low.contains("growth")) && r.id.ends_with("_B")
        }))
        .ok_or_else(|| anyhow::anyhow!(
            "Cannot find biomass reaction for species B (tried id={:?} and name heuristic)",
            biomass_b_id
        ))?;
    Ok((bio_a, bio_b))
}

/// Infer a metabolite's compartment from its ID suffix.
/// Used as a fallback when the SBML `compartment` annotation is absent.
fn get_compartment_from_id(id: &str) -> &str {
    // Bracket-form compartments (BiGG legacy / COBRA Matlab) handled first.
    if      id.ends_with("_e0") || id.ends_with("_e") || id.ends_with("[e]") { "e" }
    else if id.ends_with("_c0") || id.ends_with("_c") || id.ends_with("[c]") { "c" }
    else if id.ends_with("_p0") || id.ends_with("_p") || id.ends_with("[p]") { "p" }
    else if id.ends_with("_m")  || id.ends_with("[m]") { "m" }
    else { "c" }
}

/// Annotation-aware compartment lookup.
/// Prefers the SBML `compartment` annotation; falls back to ID-suffix heuristic.
/// Returns `true` when the metabolite is extracellular by either source.
#[inline]
fn is_extracellular_met_id(met_id: &str, annotation: Option<&str>) -> bool {
    if let Some(comp) = annotation {
        model::is_extracellular_compartment(comp)
    } else {
        get_compartment_from_id(met_id) == "e"
    }
}

type MetRxnIndex = HashMap<String, Vec<(String, f64)>>;

fn build_met_to_rxns_index(model: &MetabolicModel) -> MetRxnIndex {
    let mut idx: MetRxnIndex = HashMap::with_capacity(model.metabolites.len());
    for rxn in model.reactions.values() {
        for (met, coeff) in &rxn.metabolites {
            idx.entry(met.clone()).or_default().push((rxn.id.clone(), *coeff));
        }
    }
    idx
}

macro_rules! silent_lp {
    (maximise $obj:expr, $vars:expr) => {{
        let mut lp = $vars.maximise($obj).using(highs);
        lp.set_verbose(false); lp
    }};
    (minimise $obj:expr, $vars:expr) => {{
        let mut lp = $vars.minimise($obj).using(highs);
        lp.set_verbose(false); lp
    }};
}

// ============================================================
// Single-species: max biomass (FBA)
// ============================================================

fn solve_fba_max_biomass(
    reactions: &[Reaction], sparse_s: &SparseS, bio_idx: usize, params: &AnalysisParams,
) -> anyhow::Result<FBAResult> {
    let mut vars = ProblemVariables::new();
    let flux_vars: Vec<_> = reactions.iter().map(|r| create_flux_var(&mut vars, r, params)).collect();
    let objective = flux_vars[bio_idx];
    let mut lp = silent_lp!(maximise objective, vars);
    lp = add_mass_balance_constraints(lp, sparse_s, &flux_vars);
    let solution = lp.solve()?;
    let obj_value = solution.value(flux_vars[bio_idx]);
    let mut fd = HashMap::new();
    for (i, r) in reactions.iter().enumerate() {
        fd.insert(r.id.clone(), solution.value(flux_vars[i]));
    }
    Ok(FBAResult { objective_value: obj_value, flux_distribution: fd, fba_optimal: obj_value })
}

// ============================================================
// Combined CycleFreeFlux + pFBA
//
// `extra_locks` lets the caller pin additional reaction fluxes at
// their FBA-pass values during the parsimonious / cycle-removal LP.
// In pairwise co-culture mode this is used to pin each species'
// biomass to its phase-1 value while interior cycles are removed.
// ============================================================

fn cycle_free_pfba(
    reactions: &[Reaction], sparse_s: &SparseS, bio_idx: usize,
    fba_fluxes: &HashMap<String, f64>, params: &AnalysisParams,
    extra_locks: &[String],
) -> anyhow::Result<Option<FBAResult>> {
    let fba_biomass = fba_fluxes.get(&reactions[bio_idx].id).copied().unwrap_or(0.0);
    let lock_set: HashSet<&str> = extra_locks.iter().map(|s| s.as_str()).collect();

    let mut vars = ProblemVariables::new();
    let flux_vars: Vec<_> = reactions.iter().map(|r| create_flux_var(&mut vars, r, params)).collect();

    // Don't penalise abs flux for: exchange reactions, the main biomass,
    // or any user-specified locked reactions (their values are pinned).
    let abs_vars: Vec<Option<Variable>> = reactions.iter().enumerate()
        .map(|(j, r)| {
            if is_exchange_id(&r.id)
                || j == bio_idx
                || lock_set.contains(r.id.as_str())
            { None }
            else { Some(vars.add(variable().min(0.0))) }
        }).collect();

    let objective: Expression = abs_vars.iter()
        .filter_map(|av| av.as_ref().map(|v| 1.0 * *v)).sum();

    let mut lp = silent_lp!(minimise objective, vars);
    lp = add_mass_balance_constraints(lp, sparse_s, &flux_vars);

    // Lock all exchanges to their FBA values.
    for (j, r) in reactions.iter().enumerate() {
        if is_exchange_id(&r.id) {
            let v_star = fba_fluxes.get(&r.id).copied().unwrap_or(0.0);
            let e: Expression = 1.0 * flux_vars[j];
            lp = lp.with(e.clone().geq(v_star - params.lock_tol));
            lp = lp.with(e.leq(v_star + params.lock_tol));
        }
    }

    // Lock any caller-supplied reactions (e.g. per-species biomasses in
    // pairwise co-culture).  Skip the main biomass reaction since it has
    // its own cap/floor constraint below.
    if !lock_set.is_empty() {
        for (j, r) in reactions.iter().enumerate() {
            if j == bio_idx { continue; }
            if !lock_set.contains(r.id.as_str()) { continue; }
            let v_star = fba_fluxes.get(&r.id).copied().unwrap_or(0.0);
            let e: Expression = 1.0 * flux_vars[j];
            lp = lp.with(e.clone().geq(v_star - params.lock_tol));
            lp = lp.with(e.leq(v_star + params.lock_tol));
        }
    }

    // Constrain biomass to its FBA optimum within the flux-lock tolerance ε
    // (band form  v* − ε ≤ v_bio ≤ v* + ε), so the parsimonious / loop-removal
    // step cannot trade growth rate for lower total flux. The FBA solution
    // (biomass v*, exchanges at their optimum) is always feasible here, so this
    // never makes the LP infeasible; it only forbids a spurious biomass drop.
    let bio_cap: Expression = 1.0 * flux_vars[bio_idx];
    lp = lp.with(bio_cap.leq(fba_biomass + params.lock_tol));
    let bio_floor: Expression = 1.0 * flux_vars[bio_idx];
    lp = lp.with(bio_floor.geq(fba_biomass - params.lock_tol));

    for (j, av_opt) in abs_vars.iter().enumerate() {
        if let Some(abs_v) = av_opt {
            let c1: Expression = 1.0 * *abs_v + (-1.0) * flux_vars[j];
            lp = lp.with(c1.geq(0.0));
            let c2: Expression = 1.0 * *abs_v + 1.0 * flux_vars[j];
            lp = lp.with(c2.geq(0.0));
        }
    }

    match lp.solve() {
        Ok(sol) => {
            let new_growth = sol.value(flux_vars[bio_idx]);
            let mut fd = HashMap::new();
            for (i, r) in reactions.iter().enumerate() {
                fd.insert(r.id.clone(), sol.value(flux_vars[i]));
            }
            Ok(Some(FBAResult { objective_value: new_growth, flux_distribution: fd, fba_optimal: fba_biomass }))
        }
        Err(_) => Ok(None),
    }
}

/// Single-species or pairwise co-culture FBA + CFF + pFBA. Equivalent to
/// `run_fba_locked(model, medium, params, &[])`.
pub fn run_fba(
    model: &MetabolicModel, medium: &HashSet<String>, params: &AnalysisParams,
) -> anyhow::Result<FBAResult> {
    run_fba_locked(model, medium, params, &[])
}

/// Same as `run_fba`, but `extra_locks` lists reaction IDs whose fluxes
/// should be additionally pinned at their FBA-pass values during the
/// CFF / pFBA pass.  Used by the pairwise co-culture code to lock per-
/// species biomass values so that minimising Σ |v_i| cannot silently
/// rebalance individual μ_i.
pub fn run_fba_locked(
    model: &MetabolicModel, _medium: &HashSet<String>, params: &AnalysisParams,
    extra_locks: &[String],
) -> anyhow::Result<FBAResult> {
    let reactions = sorted_reactions(model);
    let metabolites = sorted_metabolites(model);
    let bio_idx = find_biomass_reaction(&reactions, model.biomass_reaction.as_deref())
        .ok_or_else(|| anyhow::anyhow!(
            "Cannot find biomass reaction in model '{}' (tried id={:?} and name heuristic). \
             Set model.biomass_reaction or ensure a reaction ID/name contains 'biomass' or 'growth'.",
            model.id, model.biomass_reaction
        ))?;
    let sparse_s = SparseS::build(&reactions, &metabolites);

    let fba = solve_fba_max_biomass(&reactions, &sparse_s, bio_idx, params)?;
    if fba.objective_value < MIN_VIABLE_GROWTH { return Ok(fba); }

    match cycle_free_pfba(&reactions, &sparse_s, bio_idx, &fba.flux_distribution, params, extra_locks)? {
        Some(refined) => {
            let drop = fba.objective_value - refined.objective_value;
            if drop > 1e-3 {
                eprintln!("    [CFF] biomass {:.4} → {:.4} (Δ={:+.4} cycle removed; {} extra locks)",
                    fba.objective_value, refined.objective_value, -drop, extra_locks.len());
            }
            Ok(refined)
        }
        None => Ok(fba),
    }
}

// ============================================================
// Merged model (pairwise)
// ============================================================

#[derive(Debug, Clone)]
struct MergedModel {
    reactions: Vec<Reaction>,
    metabolites: Vec<Metabolite>,
    extracellular_metabolite_ids: Vec<String>,
    biomass_a_id: Option<String>,
    biomass_b_id: Option<String>,
}

fn build_merged_model(model_a: &MetabolicModel, model_b: &MetabolicModel) -> MergedModel {
    // Detect each model's external compartment once — annotation-first, then
    // boundary-reaction count, then "e" fallback (see medium::find_external_compartment).
    let ext_comp_a = medium::find_external_compartment(model_a);
    let ext_comp_b = medium::find_external_compartment(model_b);

    let mut reactions = Vec::new();
    let mut metabolite_map: HashMap<String, Metabolite> = HashMap::new();
    let mut shared_exchange_ids = HashSet::new();

    // Insert a metabolite into the merged map.
    // `compartment_hint` is the SBML annotation from the source model; if absent
    // we fall back to ID-suffix heuristic (`get_compartment_from_id`).
    let add_met = |id: &str, suffix: &str,
                   map: &mut HashMap<String, Metabolite>,
                   compartment_hint: Option<&str>| {
        let full_id = if suffix.is_empty() {
            id.to_string()
        } else {
            format!("{}_{}", id, suffix)
        };
        map.entry(full_id.clone()).or_insert_with(|| {
            let comp = compartment_hint
                .map(|c| c.to_string())
                .unwrap_or_else(|| get_compartment_from_id(id).to_string());
            Metabolite { id: full_id, name: id.to_string(), compartment: Some(comp), boundary: false, formula: None }
        });
    };

    let mut a_ids: Vec<_> = model_a.reactions.keys().cloned().collect();
    a_ids.sort();
    for id in &a_ids {
        let r = &model_a.reactions[id];
        // Use the canonical COBRApy-style check (annotation-first) to decide
        // whether this reaction is a true environment↔cell exchange that should
        // be shared between the two species in the merged model.
        if model::is_canonical_exchange(r, model_a, &ext_comp_a) {
            if !shared_exchange_ids.contains(id) {
                let mut new_rxn = r.clone();
                new_rxn.metabolites = r.metabolites.iter().map(|(m, c)| {
                    let ann = model_a.metabolites.get(m).and_then(|mt| mt.compartment.as_deref());
                    add_met(m, "", &mut metabolite_map, ann);
                    (m.clone(), *c)
                }).collect();
                if let Some(b_rxn) = model_b.reactions.get(id) {
                    new_rxn.lower_bound = new_rxn.lower_bound.min(b_rxn.lower_bound);
                    new_rxn.upper_bound = new_rxn.upper_bound.max(b_rxn.upper_bound);
                    if b_rxn.reversible { new_rxn.reversible = true; }
                }
                reactions.push(new_rxn);
                shared_exchange_ids.insert(id.clone());
            }
        } else {
            let mut new_rxn = r.clone();
            new_rxn.id = format!("{}_A", id);
            new_rxn.metabolites = r.metabolites.iter().map(|(m, c)| {
                let ann = model_a.metabolites.get(m).and_then(|mt| mt.compartment.as_deref());
                if is_extracellular_met_id(m, ann) {
                    add_met(m, "", &mut metabolite_map, ann);
                    (m.clone(), *c)
                } else {
                    let nm = format!("{}_A", m);
                    add_met(m, "A", &mut metabolite_map, ann);
                    (nm, *c)
                }
            }).collect();
            reactions.push(new_rxn);
        }
    }

    let mut b_ids: Vec<_> = model_b.reactions.keys().cloned().collect();
    b_ids.sort();
    for id in &b_ids {
        let r = &model_b.reactions[id];
        if model::is_canonical_exchange(r, model_b, &ext_comp_b) {
            if !shared_exchange_ids.contains(id) {
                let mut new_rxn = r.clone();
                new_rxn.metabolites = r.metabolites.iter().map(|(m, c)| {
                    let ann = model_b.metabolites.get(m).and_then(|mt| mt.compartment.as_deref());
                    add_met(m, "", &mut metabolite_map, ann);
                    (m.clone(), *c)
                }).collect();
                reactions.push(new_rxn);
                shared_exchange_ids.insert(id.clone());
            }
        } else {
            let mut new_rxn = r.clone();
            new_rxn.id = format!("{}_B", id);
            new_rxn.metabolites = r.metabolites.iter().map(|(m, c)| {
                let ann = model_b.metabolites.get(m).and_then(|mt| mt.compartment.as_deref());
                if is_extracellular_met_id(m, ann) {
                    add_met(m, "", &mut metabolite_map, ann);
                    (m.clone(), *c)
                } else {
                    let nm = format!("{}_B", m);
                    add_met(m, "B", &mut metabolite_map, ann);
                    (nm, *c)
                }
            }).collect();
            reactions.push(new_rxn);
        }
    }

    // A metabolite belongs to the shared extracellular pool when its compartment
    // annotation (set by `add_met` above) is recognised as external.
    let mut ext_ids: Vec<String> = metabolite_map.values()
        .filter(|m| m.compartment.as_deref()
            .map(model::is_extracellular_compartment)
            .unwrap_or(false))
        .map(|m| m.id.clone())
        .collect();
    ext_ids.sort();
    let mut metabolites: Vec<Metabolite> = metabolite_map.into_values().collect();
    metabolites.sort_by(|a, b| a.id.cmp(&b.id));

    MergedModel {
        reactions, metabolites, extracellular_metabolite_ids: ext_ids,
        biomass_a_id: model_a.biomass_reaction.clone(),
        biomass_b_id: model_b.biomass_reaction.clone(),
    }
}

fn compute_species_target_contribution(
    merged: &MergedModel, co_fluxes: &HashMap<String, f64>,
    ext_metabolite_id: &str, species_suffix: &str,
) -> f64 {
    let mut total = 0.0;
    for rxn in &merged.reactions {
        if !rxn.id.ends_with(species_suffix) { continue; }
        for (met_id, coeff) in &rxn.metabolites {
            if met_id == ext_metabolite_id {
                let flux = co_fluxes.get(&rxn.id).copied().unwrap_or(0.0);
                total += coeff * flux;
            }
        }
    }
    total
}

fn find_exchange_metabolite(merged: &MergedModel, rxn_id: &str) -> Option<String> {
    merged.reactions.iter().find(|r| r.id == rxn_id)
        .and_then(|r| r.metabolites.first().map(|(m, _)| m.clone()))
}

// ============================================================
// Pairwise co-culture: lex max-min, then CFF+pFBA
// ============================================================

fn lp_max_min_ratio(
    merged: &MergedModel, sparse_s: &SparseS,
    growth_a_alone: f64, growth_b_alone: f64,
    a_viable: bool, b_viable: bool, params: &AnalysisParams,
) -> anyhow::Result<f64> {
    let (bio_a, bio_b) = find_merged_biomass_indices(
        &merged.reactions, merged.biomass_a_id.as_deref(), merged.biomass_b_id.as_deref())?;

    let mut vars = ProblemVariables::new();
    let flux_vars: Vec<_> = merged.reactions.iter()
        .map(|r| create_coculture_flux_var(&mut vars, r, a_viable, b_viable, params)).collect();
    let z = vars.add(variable().min(0.0));

    let mut lp = silent_lp!(maximise 1.0 * z, vars);
    lp = add_mass_balance_constraints(lp, sparse_s, &flux_vars);

    if a_viable && growth_a_alone > MIN_VIABLE_GROWTH {
        let c: Expression = growth_a_alone * z + (-1.0) * flux_vars[bio_a];
        lp = lp.with(c.leq(0.0));
    }
    if b_viable && growth_b_alone > MIN_VIABLE_GROWTH {
        let c: Expression = growth_b_alone * z + (-1.0) * flux_vars[bio_b];
        lp = lp.with(c.leq(0.0));
    }

    let sol = lp.solve()?;
    Ok(sol.value(z))
}

fn lp_max_total_at_z(
    merged: &MergedModel, sparse_s: &SparseS,
    growth_a_alone: f64, growth_b_alone: f64, z_star: f64,
    a_viable: bool, b_viable: bool, params: &AnalysisParams,
) -> anyhow::Result<(f64, f64, HashMap<String, f64>)> {
    let (bio_a, bio_b) = find_merged_biomass_indices(
        &merged.reactions, merged.biomass_a_id.as_deref(), merged.biomass_b_id.as_deref())?;

    let mut vars = ProblemVariables::new();
    let flux_vars: Vec<_> = merged.reactions.iter()
        .map(|r| create_coculture_flux_var(&mut vars, r, a_viable, b_viable, params)).collect();

    let objective: Expression = 1.0 * flux_vars[bio_a] + 1.0 * flux_vars[bio_b];
    let mut lp = silent_lp!(maximise objective, vars);
    lp = add_mass_balance_constraints(lp, sparse_s, &flux_vars);

    let z_floor = (z_star - NUMERICAL_TOL).max(0.0);
    if a_viable && growth_a_alone > MIN_VIABLE_GROWTH {
        let c: Expression = 1.0 * flux_vars[bio_a];
        lp = lp.with(c.geq(z_floor * growth_a_alone));
    }
    if b_viable && growth_b_alone > MIN_VIABLE_GROWTH {
        let c: Expression = 1.0 * flux_vars[bio_b];
        lp = lp.with(c.geq(z_floor * growth_b_alone));
    }

    let sol = lp.solve()?;
    let ga = sol.value(flux_vars[bio_a]);
    let gb = sol.value(flux_vars[bio_b]);
    let mut fluxes = HashMap::new();
    for (i, rxn) in merged.reactions.iter().enumerate() {
        fluxes.insert(rxn.id.clone(), sol.value(flux_vars[i]));
    }
    Ok((ga, gb, fluxes))
}

/// Fixed-ratio co-culture allocation (objective-sensitivity alternative to
/// the lexicographic max-min).  Maximises total community biomass subject to
/// the constraint that the co-culture growth-rate ratio equals the
/// monoculture ratio, μ_A^co/μ_B^co = μ_A^alone/μ_B^alone, applied only when
/// both partners are viable in monoculture.
fn lp_fixed_ratio(
    merged: &MergedModel, sparse_s: &SparseS,
    growth_a_alone: f64, growth_b_alone: f64,
    a_viable: bool, b_viable: bool, params: &AnalysisParams,
) -> anyhow::Result<(f64, f64, HashMap<String, f64>)> {
    let (bio_a, bio_b) = find_merged_biomass_indices(
        &merged.reactions, merged.biomass_a_id.as_deref(), merged.biomass_b_id.as_deref())?;

    let mut vars = ProblemVariables::new();
    let flux_vars: Vec<_> = merged.reactions.iter()
        .map(|r| create_coculture_flux_var(&mut vars, r, a_viable, b_viable, params)).collect();

    let objective: Expression = 1.0 * flux_vars[bio_a] + 1.0 * flux_vars[bio_b];
    let mut lp = silent_lp!(maximise objective, vars);
    lp = add_mass_balance_constraints(lp, sparse_s, &flux_vars);

    // g_b_alone * bio_a − g_a_alone * bio_b = 0  ⇔  bio_a/bio_b = g_a_alone/g_b_alone
    if a_viable && b_viable
        && growth_a_alone > MIN_VIABLE_GROWTH && growth_b_alone > MIN_VIABLE_GROWTH {
        let c: Expression =
            growth_b_alone * flux_vars[bio_a] + (-growth_a_alone) * flux_vars[bio_b];
        lp = lp.with(c.clone().leq(0.0));
        lp = lp.with(c.geq(0.0));
    }

    let sol = lp.solve()?;
    let ga = sol.value(flux_vars[bio_a]);
    let gb = sol.value(flux_vars[bio_b]);
    let mut fluxes = HashMap::new();
    for (i, rxn) in merged.reactions.iter().enumerate() {
        fluxes.insert(rxn.id.clone(), sol.value(flux_vars[i]));
    }
    Ok((ga, gb, fluxes))
}

fn cycle_free_pfba_coculture(
    merged: &MergedModel, sparse_s: &SparseS,
    fluxes: &HashMap<String, f64>,
    a_viable: bool, b_viable: bool, params: &AnalysisParams,
) -> anyhow::Result<Option<(f64, f64, HashMap<String, f64>)>> {
    let (bio_a, bio_b) = find_merged_biomass_indices(
        &merged.reactions, merged.biomass_a_id.as_deref(), merged.biomass_b_id.as_deref())?;
    let bio_a_value = fluxes.get(&merged.reactions[bio_a].id).copied().unwrap_or(0.0);
    let bio_b_value = fluxes.get(&merged.reactions[bio_b].id).copied().unwrap_or(0.0);

    let mut vars = ProblemVariables::new();
    let flux_vars: Vec<_> = merged.reactions.iter()
        .map(|r| create_coculture_flux_var(&mut vars, r, a_viable, b_viable, params)).collect();

    let abs_vars: Vec<Option<Variable>> = merged.reactions.iter().enumerate()
        .map(|(j, r)| {
            let is_exch = is_exchange_id(&r.id);
            let is_bio  = j == bio_a || j == bio_b;
            if is_exch || is_bio { None }
            else { Some(vars.add(variable().min(0.0))) }
        }).collect();

    let objective: Expression = abs_vars.iter()
        .filter_map(|av| av.as_ref().map(|v| 1.0 * *v)).sum();

    let mut lp = silent_lp!(minimise objective, vars);
    lp = add_mass_balance_constraints(lp, sparse_s, &flux_vars);

    for (j, r) in merged.reactions.iter().enumerate() {
        if is_exchange_id(&r.id) {
            let v_star = fluxes.get(&r.id).copied().unwrap_or(0.0);
            let e: Expression = 1.0 * flux_vars[j];
            lp = lp.with(e.clone().geq(v_star - params.lock_tol));
            lp = lp.with(e.leq(v_star + params.lock_tol));
        }
    }

    let bio_a_cap: Expression = 1.0 * flux_vars[bio_a];
    lp = lp.with(bio_a_cap.leq(bio_a_value));
    let bio_b_cap: Expression = 1.0 * flux_vars[bio_b];
    lp = lp.with(bio_b_cap.leq(bio_b_value));

    for (j, av_opt) in abs_vars.iter().enumerate() {
        if let Some(abs_v) = av_opt {
            let c1: Expression = 1.0 * *abs_v + (-1.0) * flux_vars[j];
            lp = lp.with(c1.geq(0.0));
            let c2: Expression = 1.0 * *abs_v + 1.0 * flux_vars[j];
            lp = lp.with(c2.geq(0.0));
        }
    }

    match lp.solve() {
        Ok(sol) => {
            let ga = sol.value(flux_vars[bio_a]);
            let gb = sol.value(flux_vars[bio_b]);
            let mut new_fluxes = HashMap::new();
            for (i, r) in merged.reactions.iter().enumerate() {
                new_fluxes.insert(r.id.clone(), sol.value(flux_vars[i]));
            }
            Ok(Some((ga, gb, new_fluxes)))
        }
        Err(_) => Ok(None),
    }
}

fn coculture_z(g_a_alone: f64, g_b_alone: f64, ga: f64, gb: f64) -> f64 {
    let r_a = if g_a_alone > MIN_VIABLE_GROWTH { ga / g_a_alone } else { f64::INFINITY };
    let r_b = if g_b_alone > MIN_VIABLE_GROWTH { gb / g_b_alone } else { f64::INFINITY };
    r_a.min(r_b)
}

fn run_co_culture(
    merged: &MergedModel, _medium: &HashSet<String>, sparse_s: &SparseS,
    growth_a_alone: f64, growth_b_alone: f64, params: &AnalysisParams,
) -> anyhow::Result<(f64, f64, f64, HashMap<String, f64>)> {
    let a_viable = growth_a_alone > MIN_VIABLE_GROWTH;
    let b_viable = growth_b_alone > MIN_VIABLE_GROWTH;

    if !a_viable && !b_viable {
        return Ok((0.0, 0.0, 0.0, HashMap::new()));
    }

    let (ga2, gb2, fluxes2) = if params.fixed_ratio {
        let r = lp_fixed_ratio(
            merged, sparse_s, growth_a_alone, growth_b_alone, a_viable, b_viable, params)?;
        eprintln!("    FixedRatio:  →  A={:.4} B={:.4}", r.0, r.1);
        r
    } else {
        let z_star = lp_max_min_ratio(
            merged, sparse_s, growth_a_alone, growth_b_alone, a_viable, b_viable, params)?;
        let r = lp_max_total_at_z(
            merged, sparse_s, growth_a_alone, growth_b_alone, z_star,
            a_viable, b_viable, params)?;
        eprintln!("    LexMaxMin: z*={:.4}  →  A={:.4} B={:.4}", z_star, r.0, r.1);
        r
    };

    let (ga, gb, fluxes) = if a_viable || b_viable {
        match cycle_free_pfba_coculture(
            merged, sparse_s, &fluxes2, a_viable, b_viable, params)?
        {
            Some((ga3, gb3, fluxes3)) => {
                let total_drop = (ga2 + gb2) - (ga3 + gb3);
                if total_drop > 1e-3 {
                    eprintln!("    [CFF] co-culture biomass {:.4} → {:.4} (Δ={:+.4})",
                        ga2 + gb2, ga3 + gb3, -total_drop);
                }
                (ga3, gb3, fluxes3)
            }
            None => (ga2, gb2, fluxes2),
        }
    } else {
        (ga2, gb2, fluxes2)
    };

    let z_post = coculture_z(growth_a_alone, growth_b_alone, ga, gb);
    Ok((ga, gb, z_post, fluxes))
}

// ============================================================
// Cross-feeding analysis
// ============================================================

fn find_genes_for_metabolite(
    ext_met_id: &str, model: &MetabolicModel, index: &MetRxnIndex, direction: f64,
) -> Vec<String> {
    let mut genes = Vec::new();
    if let Some(entries) = index.get(ext_met_id) {
        for (rxn_id, coeff) in entries {
            if *coeff * direction <= 0.0 { continue; }
            if let Some(rxn) = model.reactions.get(rxn_id) {
                for gid in &rxn.gene_ids {
                    let label = model.genes.get(gid)
                        .map(|g| if g.label.is_empty() { g.id.clone() } else { g.label.clone() })
                        .unwrap_or_else(|| gid.clone());
                    genes.push(label);
                }
            }
        }
    }
    genes.sort(); genes.dedup();
    genes
}

fn analyze_cross_feeding(
    merged: &MergedModel, sparse_s: &SparseS, fluxes: &HashMap<String, f64>,
    model_a: &MetabolicModel, model_b: &MetabolicModel,
    idx_a: &MetRxnIndex, idx_b: &MetRxnIndex,
    growth_a_co: f64, growth_b_co: f64,
    fba_a: &FBAResult, fba_b: &FBAResult,
    growth_a_alone: f64, growth_b_alone: f64, _params: &AnalysisParams,
) -> (Vec<MetaboliteExchange>, Vec<MetaboliteExchange>, Vec<SharedResource>) {
    let mut a_to_b = Vec::new();
    let mut b_to_a = Vec::new();
    let mut shared = Vec::new();

    let inorganics      = inorganic_set();
    let cofactors       = cofactor_set();
    let reactive_inters = reactive_intermediate_set();
    let electron_carrs  = electron_carrier_set();
    let low_conf_set    = low_confidence_crossfeed_set();
    // Build id -> metabolite-name lookup so we can run the name-based
    // CoA-derivative detector on each extracellular metabolite.
    let met_name: HashMap<&str, &str> = merged.metabolites.iter()
        .map(|m| (m.id.as_str(), m.name.as_str()))
        .collect();

    let mut single_a_ex: HashMap<String, f64> = HashMap::new();
    let mut single_b_ex: HashMap<String, f64> = HashMap::new();
    for rxn_id in model_a.reactions.keys() {
        if is_exchange_id(rxn_id) {
            single_a_ex.insert(rxn_id.clone(),
                fba_a.flux_distribution.get(rxn_id).copied().unwrap_or(0.0));
        }
    }
    for rxn_id in model_b.reactions.keys() {
        if is_exchange_id(rxn_id) {
            single_b_ex.insert(rxn_id.clone(),
                fba_b.flux_distribution.get(rxn_id).copied().unwrap_or(0.0));
        }
    }

    let mut direct_cf_mets: HashSet<String> = HashSet::new();

    struct MetContrib {
        a_contrib: f64, b_contrib: f64,
        is_inorganic: bool, bound_saturated: bool,
    }
    let mut met_contribs: HashMap<String, MetContrib> = HashMap::new();
    let met_row_idx: HashMap<&str, usize> = merged.metabolites.iter().enumerate()
        .map(|(i, m)| (m.id.as_str(), i)).collect();

    for ext_met_id in &merged.extracellular_metabolite_ids {
        let mut a_contrib = 0.0f64;
        let mut b_contrib = 0.0f64;
        let mut ex_contrib = 0.0f64;
        let mut bound_saturated = false;

        if let Some(&row_i) = met_row_idx.get(ext_met_id.as_str()) {
            for &(j, coeff) in &sparse_s.rows[row_i] {
                let rxn = &merged.reactions[j];
                let flux = fluxes.get(&rxn.id).copied().unwrap_or(0.0);
                if flux.abs() < NUMERICAL_TOL { continue; }
                if rxn.lower_bound < -NUMERICAL_TOL
                    && (flux - rxn.lower_bound).abs() < NUMERICAL_TOL {
                    bound_saturated = true;
                }
                if rxn.upper_bound > NUMERICAL_TOL
                    && (flux - rxn.upper_bound).abs() < NUMERICAL_TOL {
                    bound_saturated = true;
                }
                let contribution = coeff * flux;
                if rxn.id.ends_with("_A") { a_contrib += contribution; }
                else if rxn.id.ends_with("_B") { b_contrib += contribution; }
                else { ex_contrib += contribution; }
            }
        }

        FILTER_EVALUATED.fetch_add(1, Ordering::Relaxed);

        let base = metabolite_base(ext_met_id);
        let is_inorganic     = inorganics.contains(&base);
        let is_cofactor      = cofactors.contains(&base);
        let is_reactive      = reactive_inters.contains(&base);
        let is_electron      = electron_carrs.contains(&base);
        // Catch the ~900 acyl-CoA derivatives by name (-CoA suffix).
        let is_coa_deriv     = met_name.get(ext_met_id.as_str())
            .map(|n| is_coa_derivative_name(n)).unwrap_or(false);
        let is_excluded      = is_inorganic || is_cofactor || is_reactive
                             || is_electron  || is_coa_deriv;
        // Low-confidence: NOT excluded (kept in output) but flagged so
        // downstream figures can render Pyruvate / Lactate / Ethanol /
        // Glycerol differently from well-established cross-feeds.
        let is_lc            = low_conf_set.contains(&base);

        // Per-category counters (attribute to first matching category so
        // each evaluation contributes to exactly one filtered/kept bucket).
        if is_inorganic         { FILTER_INORGANIC.fetch_add(1, Ordering::Relaxed); }
        else if is_cofactor     { FILTER_COFACTOR.fetch_add(1, Ordering::Relaxed); }
        else if is_coa_deriv    { FILTER_COA_DERIV.fetch_add(1, Ordering::Relaxed); }
        else if is_reactive     { FILTER_REACTIVE.fetch_add(1, Ordering::Relaxed); }
        else if is_electron     { FILTER_ELECTRON.fetch_add(1, Ordering::Relaxed); }
        else if bound_saturated { FILTER_BOUND_SATURATED.fetch_add(1, Ordering::Relaxed); }
        else                    { FILTER_KEPT.fetch_add(1, Ordering::Relaxed); }

        met_contribs.insert(ext_met_id.clone(), MetContrib {
            a_contrib, b_contrib, is_inorganic: is_excluded, bound_saturated,
        });

        if is_excluded { continue; }
        if bound_saturated { continue; }

        if a_contrib > NUMERICAL_TOL && b_contrib < -NUMERICAL_TOL {
            let flux_val = a_contrib.min(-b_contrib);
            let dg = find_genes_for_metabolite(ext_met_id, model_a, idx_a, 1.0);
            let rg = find_genes_for_metabolite(ext_met_id, model_b, idx_b, -1.0);
            a_to_b.push(MetaboliteExchange {
                metabolite_id: ext_met_id.clone(), flux: flux_val,
                donor_genes: dg, receiver_genes: rg, is_essential: false,
                inferred: false,
                low_confidence: is_lc,
            });
            direct_cf_mets.insert(ext_met_id.clone());
        }

        if b_contrib > NUMERICAL_TOL && a_contrib < -NUMERICAL_TOL {
            let flux_val = b_contrib.min(-a_contrib);
            let dg = find_genes_for_metabolite(ext_met_id, model_b, idx_b, 1.0);
            let rg = find_genes_for_metabolite(ext_met_id, model_a, idx_a, -1.0);
            b_to_a.push(MetaboliteExchange {
                metabolite_id: ext_met_id.clone(), flux: flux_val,
                donor_genes: dg, receiver_genes: rg, is_essential: false,
                inferred: false,
                low_confidence: is_lc,
            });
            direct_cf_mets.insert(ext_met_id.clone());
        }

        let consuming_a = a_contrib < -NUMERICAL_TOL;
        let consuming_b = b_contrib < -NUMERICAL_TOL;
        if consuming_a && consuming_b {
            let a_cons = a_contrib.abs(); let b_cons = b_contrib.abs();
            let total_cons = a_cons + b_cons;
            shared.push(SharedResource {
                metabolite_id: ext_met_id.clone(),
                uptake_a: a_cons, uptake_b: b_cons, total_available: ex_contrib,
                competition_ratio: if total_cons > NUMERICAL_TOL { a_cons / total_cons } else { 0.5 },
                genes_a: find_genes_for_metabolite(ext_met_id, model_a, idx_a, -1.0),
                genes_b: find_genes_for_metabolite(ext_met_id, model_b, idx_b, -1.0),
            });
        }
    }

    let benefit_a_est = if growth_a_alone > NUMERICAL_TOL {
        (growth_a_co - growth_a_alone) / growth_a_alone
    } else { 0.0 };
    let benefit_b_est = if growth_b_alone > NUMERICAL_TOL {
        (growth_b_co - growth_b_alone) / growth_b_alone
    } else { 0.0 };

    if benefit_b_est > 1e-3 && a_to_b.is_empty() {
        let mut cands: Vec<(String, f64)> = Vec::new();
        for ext_met_id in &merged.extracellular_metabolite_ids {
            if direct_cf_mets.contains(ext_met_id) { continue; }
            let mc = match met_contribs.get(ext_met_id) { Some(m) => m, None => continue };
            if mc.is_inorganic || mc.bound_saturated { continue; }
            if mc.b_contrib >= -NUMERICAL_TOL { continue; }
            if mc.a_contrib < -NUMERICAL_TOL { continue; }

            for rxn in &merged.reactions {
                if !is_exchange_id(&rxn.id) { continue; }
                if !rxn.metabolites.iter().any(|(m,_)| m == ext_met_id) { continue; }
                let b_alone = single_b_ex.get(&rxn.id).copied().unwrap_or(0.0);
                let b_consumed_alone = b_alone < -NUMERICAL_TOL;
                let b_uptake_co = (-mc.b_contrib).abs();
                if mc.a_contrib >= 0.0 {
                    let f = if mc.a_contrib > NUMERICAL_TOL { mc.a_contrib.min(b_uptake_co) }
                        else if b_consumed_alone {
                            let inc = b_uptake_co - b_alone.abs();
                            if inc > NUMERICAL_TOL { inc } else { 0.0 }
                        } else { b_uptake_co };
                    if f > NUMERICAL_TOL { cands.push((ext_met_id.clone(), f)); }
                }
                break;
            }
        }
        cands.sort_by(|a,b| b.1.total_cmp(&a.1));
        for (m, f) in &cands {
            let dg = find_genes_for_metabolite(m, model_a, idx_a, 1.0);
            let rg = find_genes_for_metabolite(m, model_b, idx_b, -1.0);
            let is_lc = low_conf_set.contains(&metabolite_base(m));
            a_to_b.push(MetaboliteExchange {
                metabolite_id: m.clone(), flux: *f, donor_genes: dg, receiver_genes: rg,
                is_essential: false,
                inferred: true,
                low_confidence: is_lc,
            });
        }
    }

    if benefit_a_est > 1e-3 && b_to_a.is_empty() {
        let mut cands: Vec<(String, f64)> = Vec::new();
        for ext_met_id in &merged.extracellular_metabolite_ids {
            if direct_cf_mets.contains(ext_met_id) { continue; }
            let mc = match met_contribs.get(ext_met_id) { Some(m) => m, None => continue };
            if mc.is_inorganic || mc.bound_saturated { continue; }
            if mc.a_contrib >= -NUMERICAL_TOL { continue; }
            if mc.b_contrib < -NUMERICAL_TOL { continue; }

            for rxn in &merged.reactions {
                if !is_exchange_id(&rxn.id) { continue; }
                if !rxn.metabolites.iter().any(|(m,_)| m == ext_met_id) { continue; }
                let a_alone = single_a_ex.get(&rxn.id).copied().unwrap_or(0.0);
                let a_consumed_alone = a_alone < -NUMERICAL_TOL;
                let a_uptake_co = (-mc.a_contrib).abs();
                if mc.b_contrib >= 0.0 {
                    let f = if mc.b_contrib > NUMERICAL_TOL { mc.b_contrib.min(a_uptake_co) }
                        else if a_consumed_alone {
                            let inc = a_uptake_co - a_alone.abs();
                            if inc > NUMERICAL_TOL { inc } else { 0.0 }
                        } else { a_uptake_co };
                    if f > NUMERICAL_TOL { cands.push((ext_met_id.clone(), f)); }
                }
                break;
            }
        }
        cands.sort_by(|a,b| b.1.total_cmp(&a.1));
        for (m, f) in &cands {
            let dg = find_genes_for_metabolite(m, model_b, idx_b, 1.0);
            let rg = find_genes_for_metabolite(m, model_a, idx_a, -1.0);
            let is_lc = low_conf_set.contains(&metabolite_base(m));
            b_to_a.push(MetaboliteExchange {
                metabolite_id: m.clone(), flux: *f, donor_genes: dg, receiver_genes: rg,
                is_essential: false,
                inferred: true,
                low_confidence: is_lc,
            });
        }
    }

    a_to_b.sort_by(|a,b| b.flux.total_cmp(&a.flux));
    b_to_a.sort_by(|a,b| b.flux.total_cmp(&a.flux));
    shared.sort_by(|a,b| (b.uptake_a + b.uptake_b).total_cmp(&(a.uptake_a + a.uptake_b)));

    (a_to_b, b_to_a, shared)
}

// ============================================================
// Interaction classification
// ============================================================

fn classify_interaction(
    benefit_a: f64, benefit_b: f64,
    _a_to_b: &[MetaboliteExchange], _b_to_a: &[MetaboliteExchange], _shared: &[SharedResource],
) -> InteractionType {
    let eps = 1e-3;
    let a_pos = benefit_a > eps; let b_pos = benefit_b > eps;
    let a_neg = benefit_a < -eps; let b_neg = benefit_b < -eps;
    let a_neu = !a_pos && !a_neg; let b_neu = !b_pos && !b_neg;
    if a_pos && b_pos { InteractionType::Mutualism }
    else if a_neg && b_neg { InteractionType::Competition }
    else if (a_pos && b_neg) || (a_neg && b_pos) { InteractionType::Parasitism }
    else if (a_pos && b_neu) || (b_pos && a_neu) { InteractionType::Commensalism }
    else if (a_neg && b_neu) || (b_neg && a_neu) { InteractionType::Amensalism }
    else { InteractionType::Neutral }
}

// ============================================================
// Public API
// ============================================================

pub fn calculate_pairwise_interaction(
    model_a: &MetabolicModel, model_b: &MetabolicModel,
    cached_fba_a: Option<&FBAResult>, cached_fba_b: Option<&FBAResult>,
    medium: &HashSet<String>, params: &AnalysisParams, target_reactions: &[String],
) -> anyhow::Result<PairwiseResult> {
    eprintln!("  Pairwise: {} ↔ {}", model_a.id, model_b.id);

    let fba_a: Cow<FBAResult> = match cached_fba_a {
        Some(r) => Cow::Borrowed(r),
        None => Cow::Owned(run_fba(model_a, medium, params)?),
    };
    let fba_b: Cow<FBAResult> = match cached_fba_b {
        Some(r) => Cow::Borrowed(r),
        None => Cow::Owned(run_fba(model_b, medium, params)?),
    };
    let growth_a_alone = fba_a.objective_value;
    let growth_b_alone = fba_b.objective_value;
    eprintln!("    Single (FBA+CFF+pFBA): A={:.4}, B={:.4}", growth_a_alone, growth_b_alone);

    let a_viable = growth_a_alone > MIN_VIABLE_GROWTH;
    let b_viable = growth_b_alone > MIN_VIABLE_GROWTH;

    let merged = build_merged_model(model_a, model_b);
    let sparse_s = SparseS::build(&merged.reactions, &merged.metabolites);
    let idx_a = build_met_to_rxns_index(model_a);
    let idx_b = build_met_to_rxns_index(model_b);

    let (growth_a_co, growth_b_co, z_value, co_fluxes) =
        run_co_culture(&merged, medium, &sparse_s, growth_a_alone, growth_b_alone, params)?;
    eprintln!("    Co-culture: A={:.4}, B={:.4}, z={:.4}", growth_a_co, growth_b_co, z_value);

    let (mut a_to_b, mut b_to_a, shared_uptakes) = if a_viable || b_viable {
        analyze_cross_feeding(
            &merged, &sparse_s, &co_fluxes,
            model_a, model_b, &idx_a, &idx_b,
            growth_a_co, growth_b_co,
            &fba_a, &fba_b, growth_a_alone, growth_b_alone, params)
    } else { (vec![], vec![], vec![]) };

    if !a_viable { a_to_b.clear(); }
    if !b_viable { b_to_a.clear(); }

    let benefit_a = snap_zero(if growth_a_alone > 1e-10 {
        (growth_a_co - growth_a_alone) / growth_a_alone
    } else if growth_a_co > MIN_VIABLE_GROWTH { 1.0 } else { 0.0 });
    let benefit_b = snap_zero(if growth_b_alone > 1e-10 {
        (growth_b_co - growth_b_alone) / growth_b_alone
    } else if growth_b_co > MIN_VIABLE_GROWTH { 1.0 } else { 0.0 });

    let total_alone = growth_a_alone + growth_b_alone;
    let total_co = growth_a_co + growth_b_co;
    let reciprocity_index = if total_alone > 0.0 { total_co / total_alone } else { 0.0 };
    let net_interaction = benefit_a + benefit_b;

    let total_cf_flux: f64 = a_to_b.iter().map(|e| e.flux).sum::<f64>()
        + b_to_a.iter().map(|e| e.flux).sum::<f64>();
    let cf_denominator = total_co.max(total_alone);
    let cross_feeding_score = if cf_denominator > MIN_VIABLE_GROWTH {
        (total_cf_flux / cf_denominator).max(0.0)
    } else { 0.0 };
    let gene_supported_fraction =
        compute_gene_supported_fraction(cross_feeding_score, &a_to_b, &b_to_a);
    let interaction_type = classify_interaction(benefit_a, benefit_b, &a_to_b, &b_to_a, &shared_uptakes);

    let target_fluxes: Vec<TargetFluxRecord> = target_reactions.iter().map(|rxn_id| {
        let alone_a = fba_a.flux_distribution.get(rxn_id).copied().unwrap_or(f64::NAN);
        let alone_b = fba_b.flux_distribution.get(rxn_id).copied().unwrap_or(f64::NAN);
        let co_total = co_fluxes.get(rxn_id).copied().unwrap_or(f64::NAN);
        let (co_a, co_b) = match find_exchange_metabolite(&merged, rxn_id) {
            Some(met_id) => (
                compute_species_target_contribution(&merged, &co_fluxes, &met_id, "_A"),
                compute_species_target_contribution(&merged, &co_fluxes, &met_id, "_B"),
            ),
            None => (f64::NAN, f64::NAN),
        };
        TargetFluxRecord {
            reaction_id: rxn_id.clone(), alone_a, alone_b, co_total,
            co_a_contribution: co_a, co_b_contribution: co_b,
        }
    }).collect();

    let n_shared = shared_uptakes.len();
    let n_a_only = b_to_a.len();
    let n_b_only = a_to_b.len();
    let union = n_shared + n_a_only + n_b_only;
    let niche_overlap_jaccard = if union > 0 { Some(n_shared as f64 / union as f64) } else { None };

    let mut sum_min_union = 0.0;
    let mut sum_max_union = 0.0;
    let mut total_a_uptake = 0.0;
    let mut total_b_uptake = 0.0;

    for sh in &shared_uptakes {
        sum_min_union += sh.uptake_a.min(sh.uptake_b);
        sum_max_union += sh.uptake_a.max(sh.uptake_b);
        total_a_uptake += sh.uptake_a;
        total_b_uptake += sh.uptake_b;
    }
    for ex in &b_to_a {
        sum_max_union += ex.flux.abs();
        total_a_uptake += ex.flux.abs();
    }
    for ex in &a_to_b {
        sum_max_union += ex.flux.abs();
        total_b_uptake += ex.flux.abs();
    }

    let denom_min = total_a_uptake.min(total_b_uptake);
    let niche_overlap_weighted = if denom_min > 1e-9 {
        Some((sum_min_union / denom_min).min(1.0))
    } else {
        None
    };

    let competition_intensity = if sum_max_union > 1e-9 {
        (sum_min_union / sum_max_union).clamp(0.0, 1.0)
    } else {
        0.0
    };

    Ok(PairwiseResult {
        species_a: model_a.id.clone(), species_b: model_b.id.clone(),
        growth_a_alone, growth_b_alone, growth_a_co, growth_b_co,
        mu: z_value, reciprocity_index,
        benefit_a, benefit_b, net_interaction,
        interaction_type, cross_feeding_score, gene_supported_fraction,
        a_to_b_exchanges: a_to_b, b_to_a_exchanges: b_to_a, shared_uptakes,
        target_fluxes, niche_overlap_jaccard, niche_overlap_weighted, competition_intensity,
    })
}
// =====================================================================
// Unit tests
// =====================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Gene;
    use std::collections::HashMap;

    // -----------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------
    fn mock_reaction(id: &str, mets: Vec<(&str, f64)>) -> Reaction {
        Reaction {
            id: id.to_string(),
            name: id.to_string(),
            metabolites: mets.into_iter().map(|(m, c)| (m.to_string(), c)).collect(),
            reversible: false,
            lower_bound: -1000.0,
            upper_bound: 1000.0,
            gene_reaction_rule: "".to_string(),
            gene_ids: vec![],
        }
    }

    fn mock_model() -> MetabolicModel {
        let mut m = MetabolicModel {
            id: "test_model".to_string(),
            name: None,
            metabolites: HashMap::new(),
            reactions: HashMap::new(),
            biomass_reaction: Some("bio".to_string()),
            genes: HashMap::new(),
        };
        m.metabolites.insert(
            "glc__D_e".to_string(),
            Metabolite {
                id: "glc__D_e".to_string(),
                name: "glucose".to_string(),
                compartment: Some("e".to_string()),
                boundary: false,
                formula: None,
            },
        );
        m.metabolites.insert(
            "ac_e".to_string(),
            Metabolite {
                id: "ac_e".to_string(),
                name: "acetate".to_string(),
                compartment: Some("e".to_string()),
                boundary: false,
                formula: None,
            },
        );
        m.reactions.insert(
            "bio".to_string(),
            Reaction {
                id: "bio".to_string(),
                name: "biomass".to_string(),
                metabolites: vec![],
                reversible: false,
                lower_bound: 0.0,
                upper_bound: 1000.0,
                gene_reaction_rule: "".to_string(),
                gene_ids: vec![],
            },
        );
        m.reactions.insert(
            "EX_glc__D_e".to_string(),
            mock_reaction("EX_glc__D_e", vec![("glc__D_e", -1.0)]),
        );
        m.reactions.insert(
            "EX_ac_e".to_string(),
            mock_reaction("EX_ac_e", vec![("ac_e", -1.0)]),
        );
        m
    }

    // -----------------------------------------------------------------
    // Pure logic: string / classification helpers
    // -----------------------------------------------------------------
    #[test]
    fn test_snap_zero() {
        assert_eq!(snap_zero(1e-10), 0.0);
        assert_eq!(snap_zero(-1e-10), 0.0);
        assert_eq!(snap_zero(1e-8), 1e-8); // above 1e-9 threshold
        assert_eq!(snap_zero(1e-4), 1e-4);
        assert_eq!(snap_zero(-1e-4), -1e-4);
        assert_eq!(snap_zero(0.0), 0.0);
    }

    #[test]
    fn test_is_exchange_id() {
        assert!(is_exchange_id("R_EX_glc__D_e"));
        assert!(is_exchange_id("EX_glc__D_e"));
        assert!(!is_exchange_id("R_BIOMASS"));
        assert!(!is_exchange_id("PFK"));
        assert!(!is_exchange_id(""));
    }

    #[test]
    fn test_strip_compartment_suffix() {
        assert_eq!(strip_compartment_suffix("glc__D_e"), "glc__D");
        assert_eq!(strip_compartment_suffix("glc__D_e0"), "glc__D");
        assert_eq!(strip_compartment_suffix("atp_c"), "atp");
        assert_eq!(strip_compartment_suffix("atp_c0"), "atp");
        assert_eq!(strip_compartment_suffix("h_p"), "h");
        assert_eq!(strip_compartment_suffix("h_p0"), "h");
        assert_eq!(strip_compartment_suffix("h_m"), "h");
        assert_eq!(strip_compartment_suffix("atp"), "atp");
    }

    #[test]
    fn test_metabolite_base() {
        assert_eq!(metabolite_base("M_glc__D_e"), "glc__d");
        assert_eq!(metabolite_base("glc__D_e"), "glc__d");
        assert_eq!(metabolite_base("glc__D"), "glc__d");
        assert_eq!(metabolite_base("M_ATP_c"), "atp");
    }

    #[test]
    fn test_get_compartment_from_id() {
        assert_eq!(get_compartment_from_id("glc__D_e"), "e");
        assert_eq!(get_compartment_from_id("atp_c"), "c");
        assert_eq!(get_compartment_from_id("h_p"), "p");
        assert_eq!(get_compartment_from_id("h_m"), "m");
        assert_eq!(get_compartment_from_id("atp"), "c"); // fallback default
        // BiGG legacy / COBRA Matlab bracket form.
        assert_eq!(get_compartment_from_id("glc__D[e]"), "e");
        assert_eq!(get_compartment_from_id("atp[c]"), "c");
        assert_eq!(get_compartment_from_id("h[p]"), "p");
        assert_eq!(get_compartment_from_id("nadh[m]"), "m");
    }

    #[test]
    fn test_inorganic_set() {
        let set = inorganic_set();
        // BiGG-style
        assert!(set.contains("h2o"));
        assert!(set.contains("h"));
        assert!(set.contains("na1"));
        assert!(set.contains("co2"));
        // ModelSEED cpd IDs (the gap that previously let inorganics leak)
        assert!(set.contains("cpd00001"));   // H2O
        assert!(set.contains("cpd00067"));   // H+
        assert!(set.contains("cpd00009"));   // Pi
        assert!(set.contains("cpd00013"));   // NH4
        assert!(set.contains("cpd00011"));   // CO2
        assert!(set.contains("cpd00007"));   // O2
        // Real metabolites must still pass
        assert!(!set.contains("glc__d"));
        assert!(!set.contains("cpd00029"));  // Acetate
        assert!(set.contains("cpd00012"));   // PPi: now filtered (empirical review)
    }

    #[test]
    fn test_fba_artifact_set() {
        let set = fba_artifact_set();
        // Category 1: cofactors / energy carriers
        assert!(set.contains("cpd00010"));   // CoA
        assert!(set.contains("cpd00002"));   // ATP
        assert!(set.contains("cpd00003"));   // NAD
        assert!(set.contains("cpd00005"));   // NADPH
        assert!(set.contains("cpd00015"));   // FAD
        assert!(set.contains("cpd00087"));   // THF
        assert!(set.contains("cpd00345"));   // 5mthf
        // Category 2: reactive intermediates
        assert!(set.contains("cpd00071"));   // Acetaldehyde
        assert!(set.contains("cpd00055"));   // Formaldehyde
        assert!(set.contains("cpd00448"));   // Glyceraldehyde
        assert!(set.contains("cpd00428"));   // Methylglyoxal
        assert!(set.contains("cpd00145"));   // Hydroxypyruvate
        // Category 3: electron carriers
        assert!(set.contains("cpd11451"));   // Menaquinol
        assert!(set.contains("cpd11620"));   // Reduced ferredoxin
        // BiGG aliases
        assert!(set.contains("coa"));
        assert!(set.contains("acald"));
        assert!(set.contains("fdxr"));
        // Category 4 (real cross-feeding) must NOT be filtered
        assert!(!set.contains("cpd00029"));  // Acetate
        assert!(!set.contains("cpd00141"));  // Propionate
        assert!(!set.contains("cpd00047"));  // Formate
        assert!(!set.contains("cpd00036"));  // Succinate
        assert!(!set.contains("cpd00100"));  // Glycerol  (real signal, not propionyl-CoA)
        assert!(!set.contains("cpd00154"));  // D-Xylose  (real signal, not butyryl-CoA)
        assert!(!set.contains("cpd00314"));  // D-Mannitol (not formaldehyde!)
        // PPi is filtered by the inorganic_set (added after empirical review),
        // NOT the artifact set -- test_inorganic_set already covers it.
        assert!(!set.contains("cpd00040"));  // Glyoxylate (kept by design)
        assert!(!set.contains("glx"));       // Glyoxylate BiGG alias (kept)
    }

    #[test]
    fn test_is_coa_derivative_name() {
        // Acyl-CoA derivatives — must be filtered by name
        assert!(is_coa_derivative_name("Acetyl-CoA"));
        assert!(is_coa_derivative_name("Propionyl-CoA"));
        assert!(is_coa_derivative_name("Succinyl-CoA"));
        assert!(is_coa_derivative_name("Malonyl-CoA"));
        assert!(is_coa_derivative_name("3-Hydroxybutyryl-CoA"));
        assert!(is_coa_derivative_name("Coenzyme A"));
        // Note: bare "CoA" is filtered by ID (cpd00010 / "coa" in HashSet),
        // not by this name pattern -- so we don't assert on it here.

        // Real metabolites — must NOT match
        assert!(!is_coa_derivative_name("Acetate"));
        assert!(!is_coa_derivative_name("Propionate"));
        assert!(!is_coa_derivative_name("Coenzyme Q"));
        assert!(!is_coa_derivative_name("Cocaine"));   // contains "co" but not as -CoA suffix

        // ── Whitelisted CoA-precursors / Vitamin B5 (must be kept!) ──
        assert!(!is_coa_derivative_name("Pantothenate"));
        assert!(!is_coa_derivative_name("Pantothenic acid"));
        assert!(!is_coa_derivative_name("Vitamin B5"));
        assert!(!is_coa_derivative_name("4'-Phosphopantetheine"));
        assert!(!is_coa_derivative_name("Pantetheine"));
    }

    #[test]
    fn test_low_confidence_crossfeed_set() {
        let set = low_confidence_crossfeed_set();
        // Borderline cross-feeds — kept in output but flagged
        assert!(set.contains("cpd00020"));   // Pyruvate
        assert!(set.contains("cpd00159"));   // L-Lactate
        assert!(set.contains("cpd00221"));   // D-Lactate
        assert!(set.contains("cpd00363"));   // Ethanol
        assert!(set.contains("cpd00100"));   // Glycerol
        // BiGG aliases
        assert!(set.contains("pyr"));
        assert!(set.contains("etoh"));
        // Well-established cross-feeds must NOT be flagged as low-confidence
        assert!(!set.contains("cpd00029"));  // Acetate
        assert!(!set.contains("cpd00141"));  // Propionate
        assert!(!set.contains("cpd00047"));  // Formate
        assert!(!set.contains("cpd00036"));  // Succinate
        // Artifacts (already filtered) must not appear here either
        assert!(!set.contains("cpd00010"));  // CoA
        assert!(!set.contains("cpd00001"));  // H2O
    }

    #[test]
    fn test_compartment_suffix_brackets() {
        // BiGG bracket-style compartments must also be stripped.
        assert_eq!(strip_compartment_suffix("coa[e]"), "coa");
        assert_eq!(strip_compartment_suffix("acald[c]"), "acald");
        assert_eq!(metabolite_base("M_coa[e]"), "coa");
        assert_eq!(metabolite_base("Glc__D[e]"), "glc__d");
    }

    #[test]
    fn test_has_real_gene_support() {
        assert!(has_real_gene_support(&vec!["gene1".to_string(), "gene2".to_string()]));
        assert!(!has_real_gene_support(&vec![]));
        assert!(!has_real_gene_support(&vec!["spontaneous".to_string()]));
        // "hypothetical" is NOT filtered by the current implementation
        assert!(has_real_gene_support(&vec!["unknown".to_string(), "hypothetical".to_string()]));
        assert!(has_real_gene_support(&vec!["unknown".to_string(), "real".to_string()]));
    }

    #[test]
    fn test_compute_gene_supported_fraction() {
        let exchanges = vec![
            MetaboliteExchange {
                metabolite_id: "ac_e".to_string(),
                flux: 5.0,
                donor_genes: vec!["acs".to_string()],
                receiver_genes: vec!["ackA".to_string()],
                is_essential: false,
                inferred: false,
                low_confidence: false,
            },
            MetaboliteExchange {
                metabolite_id: "lac__D_e".to_string(),
                flux: 3.0,
                donor_genes: vec![],
                receiver_genes: vec![],
                is_essential: false,
                inferred: false,
                low_confidence: true,  // D-Lactate flagged as low-confidence
            },
        ];
        // 5.0 with genes + 3.0 without = 5/8 = 0.625
        let frac = compute_gene_supported_fraction(0.5, &exchanges, &vec![]);
        assert!((frac - 0.625).abs() < 1e-6, "frac={}", frac);
    }

    #[test]
    fn test_classify_interaction() {
        let empty_ex: Vec<MetaboliteExchange> = vec![];
        let empty_sh: Vec<SharedResource> = vec![];
        assert_eq!(classify_interaction(0.1, 0.1, &empty_ex, &empty_ex, &empty_sh), InteractionType::Mutualism);
        assert_eq!(classify_interaction(0.1, -0.1, &empty_ex, &empty_ex, &empty_sh), InteractionType::Parasitism);
        assert_eq!(classify_interaction(-0.1, 0.1, &empty_ex, &empty_ex, &empty_sh), InteractionType::Parasitism);
        assert_eq!(classify_interaction(-0.1, -0.1, &empty_ex, &empty_ex, &empty_sh), InteractionType::Competition);
        assert_eq!(classify_interaction(0.1, 0.0, &empty_ex, &empty_ex, &empty_sh), InteractionType::Commensalism);
        assert_eq!(classify_interaction(0.0, 0.1, &empty_ex, &empty_ex, &empty_sh), InteractionType::Commensalism);
        assert_eq!(classify_interaction(-0.1, 0.0, &empty_ex, &empty_ex, &empty_sh), InteractionType::Amensalism);
        assert_eq!(classify_interaction(0.0, -0.1, &empty_ex, &empty_ex, &empty_sh), InteractionType::Amensalism);
        assert_eq!(classify_interaction(0.0, 0.0, &empty_ex, &empty_ex, &empty_sh), InteractionType::Neutral);
    }

    #[test]
    fn test_coculture_z() {
        assert_eq!(coculture_z(1.0, 1.0, 0.5, 0.5), 0.5);
        assert_eq!(coculture_z(1.0, 1.0, 0.0, 0.5), 0.0);
        // When g_a_alone <= MIN_VIABLE_GROWTH, r_a = INF, so min = r_b
        assert_eq!(coculture_z(0.0, 1.0, 0.5, 0.5), 0.5);
        assert_eq!(coculture_z(1.0, 2.0, 0.5, 1.0), 0.5);
    }

    // -----------------------------------------------------------------
    // Medium-difficulty: in-memory data structures
    // -----------------------------------------------------------------
    #[test]
    fn test_sorted_reactions_and_metabolites() {
        let model = mock_model();
        let rxns = sorted_reactions(&model);
        assert_eq!(rxns.len(), 3);
        assert_eq!(rxns[0].id, "EX_ac_e");
        assert_eq!(rxns[1].id, "EX_glc__D_e");
        assert_eq!(rxns[2].id, "bio");

        let mets = sorted_metabolites(&model);
        assert_eq!(mets.len(), 2);
        assert_eq!(mets[0].id, "ac_e");
        assert_eq!(mets[1].id, "glc__D_e");
    }

    #[test]
    fn test_find_biomass_reaction() {
        let rxns = vec![
            mock_reaction("PFK", vec![]),
            mock_reaction("bio", vec![]),
        ];
        assert_eq!(find_biomass_reaction(&rxns, Some("bio")), Some(1));
        // When no biomass id given and no name match, returns None (no silent fallback).
        assert_eq!(find_biomass_reaction(&rxns, None), None);

        let rxns2 = vec![
            mock_reaction("PFK", vec![]),
            mock_reaction("BIOMASS", vec![]),
        ];
        assert_eq!(find_biomass_reaction(&rxns2, None), Some(1));
    }

    #[test]
    fn test_find_merged_biomass_indices() {
        let rxns = vec![
            mock_reaction("biomass_A", vec![]),
            mock_reaction("PFK_A", vec![]),
            mock_reaction("biomass_B", vec![]),
            mock_reaction("PFK_B", vec![]),
        ];
        let (a, b) = find_merged_biomass_indices(&rxns, Some("biomass_A"), Some("biomass_B"))
            .expect("should find both biomass reactions");
        assert_eq!(a, 0);
        assert_eq!(b, 2);
    }

    #[test]
    fn test_build_met_to_rxns_index() {
        let model = mock_model();
        let index = build_met_to_rxns_index(&model);
        assert!(index.contains_key("glc__D_e"));
        assert!(index.contains_key("ac_e"));
        let glc_entries = index.get("glc__D_e").unwrap();
        assert_eq!(glc_entries.len(), 1);
        assert_eq!(glc_entries[0].0, "EX_glc__D_e");
        assert_eq!(glc_entries[0].1, -1.0);
    }

    #[test]
    fn test_build_merged_model() {
        let m1 = mock_model();
        let mut m2 = mock_model();
        m2.id = "model2".to_string();
        let merged = build_merged_model(&m1, &m2);

        // Should have shared extracellular metabolites
        assert!(merged.metabolites.iter().any(|m| m.id == "glc__D_e"));
        assert!(merged.metabolites.iter().any(|m| m.id == "ac_e"));

        // Internal reactions get _A / _B suffix
        assert!(merged.reactions.iter().any(|r| r.id == "bio_A"));
        assert!(merged.reactions.iter().any(|r| r.id == "bio_B"));

        // Exchange reactions are shared (no suffix) because both models have the same exchange
        assert!(merged.reactions.iter().any(|r| r.id == "EX_glc__D_e"));
        assert!(!merged.reactions.iter().any(|r| r.id == "EX_glc__D_e_A"));
        assert!(!merged.reactions.iter().any(|r| r.id == "EX_glc__D_e_B"));
    }

    #[test]
    fn test_find_exchange_metabolite() {
        let m1 = mock_model();
        let m2 = mock_model();
        let merged = build_merged_model(&m1, &m2);
        // Exchange reactions are shared (no suffix) in merged model
        assert_eq!(find_exchange_metabolite(&merged, "EX_glc__D_e"), Some("glc__D_e".to_string()));
        assert_eq!(find_exchange_metabolite(&merged, "bio_A"), None); // not exchange, 0 mets
        assert_eq!(find_exchange_metabolite(&merged, "nonexistent"), None);
    }

    #[test]
    fn test_compute_species_target_contribution() {
        // Build a minimal merged model with species-specific reactions containing the target metabolite
        let merged = MergedModel {
            reactions: vec![
                Reaction {
                    id: "GLCpts_A".to_string(),
                    name: "glucose transport A".to_string(),
                    metabolites: vec![("glc__D_e".to_string(), -1.0), ("glc__D_c_A".to_string(), 1.0)],
                    reversible: false,
                    lower_bound: 0.0,
                    upper_bound: 1000.0,
                    gene_reaction_rule: "".to_string(),
                    gene_ids: vec![],
                },
                Reaction {
                    id: "GLCpts_B".to_string(),
                    name: "glucose transport B".to_string(),
                    metabolites: vec![("glc__D_e".to_string(), -1.0), ("glc__D_c_B".to_string(), 1.0)],
                    reversible: false,
                    lower_bound: 0.0,
                    upper_bound: 1000.0,
                    gene_reaction_rule: "".to_string(),
                    gene_ids: vec![],
                },
            ],
            metabolites: vec![
                Metabolite { id: "glc__D_e".to_string(), name: "glucose".to_string(), compartment: Some("e".to_string()), boundary: false, formula: None },
            ],
            extracellular_metabolite_ids: vec!["glc__D_e".to_string()],
            biomass_a_id: None,
            biomass_b_id: None,
        };
        let mut fluxes = HashMap::new();
        fluxes.insert("GLCpts_A".to_string(), -5.0);
        fluxes.insert("GLCpts_B".to_string(), -3.0);

        let a = compute_species_target_contribution(&merged, &fluxes, "glc__D_e", "_A");
        let b = compute_species_target_contribution(&merged, &fluxes, "glc__D_e", "_B");
        // coeff * flux = (-1.0) * (-5.0) = 5.0
        assert_eq!(a, 5.0);
        assert_eq!(b, 3.0);
    }

    #[test]
    fn test_find_genes_for_metabolite() {
        let mut model = mock_model();
        model.reactions.insert(
            "ACS".to_string(),
            Reaction {
                id: "ACS".to_string(),
                name: "acs".to_string(),
                metabolites: vec![("ac_e".to_string(), -1.0)],
                reversible: false,
                lower_bound: 0.0,
                upper_bound: 1000.0,
                gene_reaction_rule: "".to_string(),
                gene_ids: vec!["gene_acs".to_string()],
            },
        );
        model.genes.insert(
            "gene_acs".to_string(),
            Gene {
                id: "gene_acs".to_string(),
                name: "acs".to_string(),
                label: "acsA".to_string(),
            },
        );
        let index = build_met_to_rxns_index(&model);
        let genes = find_genes_for_metabolite("ac_e", &model, &index, -1.0);
        assert!(genes.contains(&"acsA".to_string()));

        // Direction check: asking for production (1.0) when reaction consumes (-1.0 coeff) should return empty
        let genes_prod = find_genes_for_metabolite("ac_e", &model, &index, 1.0);
        assert!(genes_prod.is_empty());
    }
}

