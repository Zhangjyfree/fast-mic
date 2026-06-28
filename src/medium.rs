use crate::model::{self, MetabolicModel};
use std::collections::{HashMap, HashSet};

// ─────────────────────────────────────────────────────────────────────
// BiGG → SEED compound mapping (for gapseq model support)
// ─────────────────────────────────────────────────────────────────────

/// Load a BiGG-abbreviation → **all** SEED-cpd-IDs map from ModelSEED's compounds.tsv.
///
/// Many BiGG compounds have multiple SEED aliases (e.g. `nh4` maps to both
/// `cpd00013` (NH3) and `cpd19013` (NH4+)).  gapseq models commonly use a
/// different SEED ID than the first one listed, so we collect all aliases and
/// add every matching SEED ID to the medium.
///
/// Returns an empty map (graceful no-op) if the file is missing or malformed,
/// so CarveMe / BiGG models continue to work unchanged.
pub fn load_bigg_to_seed(path: &str) -> HashMap<String, Vec<String>> {
    let content = match std::fs::read_to_string(path) {
        Ok(c) => c,
        Err(_) => return HashMap::new(),
    };
    let mut map: HashMap<String, Vec<String>> = HashMap::new();
    let mut lines = content.lines();
    let header = lines.next().unwrap_or("");
    let cols: Vec<&str> = header.split('\t').collect();
    let id_col    = cols.iter().position(|&c| c == "id").unwrap_or(0);
    let alias_col = cols.iter().position(|&c| c == "aliases").unwrap_or(18);
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= alias_col { continue; }
        let seed_id = fields[id_col].trim();
        let aliases = fields[alias_col];
        if !aliases.contains("BiGG:") { continue; }
        for seg in aliases.split('|') {
            let seg = seg.trim();
            if let Some(rest) = seg.strip_prefix("BiGG:") {
                for b in rest.split(';') {
                    let b = b.trim();
                    if !b.is_empty() {
                        let entry = map.entry(b.to_string()).or_default();
                        if !entry.contains(&seed_id.to_string()) {
                            entry.push(seed_id.to_string());
                        }
                    }
                }
            }
        }
    }
    map
}

// ─────────────────────────────────────────────────────────────────────
// External compartment detection
// Equivalent to COBRApy's `find_external_compartment` in boundary_types.py
// ─────────────────────────────────────────────────────────────────────

/// Detect the external compartment from model metabolite annotations.
/// Mimics COBRApy `find_external_compartment(model)`.
///
/// Uses `model::KNOWN_EXTERNAL_COMPARTMENTS` as the single source of truth
/// for recognised external compartment names.
pub fn find_external_compartment(model: &MetabolicModel) -> String {
    // Collect unique compartment IDs from metabolites
    let mut compartments: HashSet<String> = HashSet::new();
    for met in model.metabolites.values() {
        if let Some(ref comp) = met.compartment {
            compartments.insert(comp.clone());
        }
    }

    // Strategy 1 (COBRApy primary): match known external names
    for comp in &compartments {
        if model::is_extracellular_compartment(comp) {
            return comp.clone();
        }
    }

    // Strategy 2 (COBRApy fallback): compartment with most boundary reactions
    let mut boundary_counts: HashMap<String, usize> = HashMap::new();
    for rxn in model.reactions.values() {
        if rxn.metabolites.len() == 1 {
            let met_id = &rxn.metabolites[0].0;
            if let Some(met) = model.metabolites.get(met_id) {
                if let Some(ref comp) = met.compartment {
                    *boundary_counts.entry(comp.clone()).or_insert(0) += 1;
                }
            }
        }
    }

    if let Some((comp, count)) = boundary_counts.iter().max_by_key(|(_, &c)| c) {
        eprintln!(
            "  Medium: external compartment detected by boundary-reaction count: '{}' ({} boundary rxns)",
            comp, count
        );
        return comp.clone();
    }

    eprintln!("  Medium: WARNING – could not detect external compartment, falling back to 'e'");
    "e".to_string()
}

// ─────────────────────────────────────────────────────────────────────
// Exchange reaction identification
// Equivalent to COBRApy's `find_boundary_types(model, "exchange")`
// Uses model::is_canonical_exchange as the single implementation.
// ─────────────────────────────────────────────────────────────────────

/// Find all exchange reactions in the model.
/// Returns a sorted Vec of reaction IDs for deterministic behaviour.
pub fn find_exchange_reactions(model: &MetabolicModel) -> Vec<String> {
    let ext_comp = find_external_compartment(model);
    let mut exchanges: Vec<String> = model
        .reactions
        .iter()
        .filter(|(_, rxn)| model::is_canonical_exchange(rxn, model, &ext_comp))
        .map(|(rid, _)| rid.clone())
        .collect();
    exchanges.sort();
    exchanges
}

// ─────────────────────────────────────────────────────────────────────
// Medium compound matching
// ─────────────────────────────────────────────────────────────────────

/// Expand medium compound names to match various SBML naming conventions.
///
/// **KEY FIX**: Only generates EXTERNAL compartment variants (`_e`, `_e0`).
/// The old version generated `_c` / `_p` variants, which could cause false
/// positive matches with cytoplasmic or periplasmic metabolites.
///
/// Given input "glc__D_e", generates:
///   glc__D_e, M_glc__D_e, glc__D, M_glc__D, glc__D_e0, M_glc__D_e0
pub fn expand_medium_compounds(base_compounds: &HashSet<String>) -> HashSet<String> {
    let mut expanded = HashSet::new();
    for raw in base_compounds {
        expanded.insert(raw.clone());

        // Strip M_ prefix if present
        let without_prefix = raw.strip_prefix("M_").unwrap_or(raw);
        expanded.insert(without_prefix.to_string());
        expanded.insert(format!("M_{}", without_prefix));

        // Strip compartment suffix to get the root compound name
        let root = without_prefix
            .strip_suffix("_e")
            .or_else(|| without_prefix.strip_suffix("_e0"))
            .or_else(|| without_prefix.strip_suffix("_c"))
            .or_else(|| without_prefix.strip_suffix("_c0"))
            .or_else(|| without_prefix.strip_suffix("_p"))
            .or_else(|| without_prefix.strip_suffix("_p0"))
            .unwrap_or(without_prefix);

        expanded.insert(root.to_string());
        expanded.insert(format!("M_{}", root));

        // Only external-compartment variants — medium compounds live outside
        for ext_suffix in &["_e", "_e0"] {
            expanded.insert(format!("{}{}", root, ext_suffix));
            expanded.insert(format!("M_{}{}", root, ext_suffix));
        }
    }
    expanded
}

/// Check whether a metabolite ID matches any expanded medium compound.
fn metabolite_in_medium(met_id: &str, expanded: &HashSet<String>) -> bool {
    // Direct match
    if expanded.contains(met_id) {
        return true;
    }
    // Without M_ prefix
    let stripped = met_id.strip_prefix("M_").unwrap_or(met_id);
    if expanded.contains(stripped) {
        return true;
    }
    // Root form (strip external-compartment suffix)
    let root = stripped
        .strip_suffix("_e")
        .or_else(|| stripped.strip_suffix("_e0"))
        .unwrap_or("");
    if !root.is_empty() {
        return expanded.contains(root) || expanded.contains(&format!("M_{}", root));
    }
    false
}

// ─────────────────────────────────────────────────────────────────────
// Compound-type classification for tiered uptake limits
// ─────────────────────────────────────────────────────────────────────

/// Tiered uptake bounds for different metabolite classes.
#[derive(Debug, Clone, Copy)]
pub struct MediumBounds {
    /// Carbon sources: sugars, organic acids, SCFAs (mmol/gDW/h)
    pub carbon_source: f64,
    /// Amino acids (mmol/gDW/h)
    pub amino_acid: f64,
    /// Nucleobases and nucleosides (mmol/gDW/h)
    pub nucleobase: f64,
    /// Vitamins and enzyme cofactors (mmol/gDW/h)
    pub cofactor: f64,
}

impl Default for MediumBounds {
    fn default() -> Self {
        Self { carbon_source: 10.0, amino_acid: 1.0, nucleobase: 0.5, cofactor: 0.1 }
    }
}

/// Return true if the metabolite is a vitamin or enzyme cofactor.
/// These are always preserved at their SBML-defined lower bound even when
/// the compound is absent from the explicit medium definition, because they
/// represent trace micronutrients that are present in every realistic growth
/// medium (vitamins, quinones, metal cofactors, coenzymes).
fn is_cofactor_metabolite(met_id: &str) -> bool {
    let lower = met_id.to_lowercase();
    let clean = lower.strip_prefix("m_").unwrap_or(&lower);
    let base = clean
        .trim_end_matches("_e").trim_end_matches("_e0")
        .trim_end_matches("_c").trim_end_matches("_p");
    const COFACTORS: &[&str] = &[
        "fol", "thf", "5mthf", "mlthf",
        "nmn", "nad", "nadp", "fmn", "fad", "ribflv",
        "cbl1", "cbl2", "adocbl", "mecbl",
        "btn", "pnto__r", "pnto__s", "4ppcys", "pap",
        "pyam5p", "pydx", "pydxn", "pydx5p",
        "thm", "thmpp",
        "pheme", "sheme", "hem",
        "mqn7", "mqn8", "mql8", "q8", "q8h2",
        "lipoate", "lipoamide",
        "spmd", "ptrc", "cadav",
    ];
    const COFACTORS_SEED: &[&str] = &[
        "cpd00003", "cpd00006", "cpd00015", "cpd00016", "cpd00028",
        "cpd00045", "cpd00050", "cpd00056", "cpd00087", "cpd00104",
        "cpd00118", "cpd00125", "cpd00166", "cpd00215", "cpd00220",
        "cpd00263", "cpd00264", "cpd00305", "cpd00355", "cpd00393",
        "cpd00423", "cpd00493", "cpd00541", "cpd00557", "cpd00635",
        "cpd00644", "cpd02666", "cpd11606", "cpd15499", "cpd15500",
        "cpd15561",
    ];
    COFACTORS.contains(&base) || COFACTORS_SEED.contains(&base)
}

/// Classify a metabolite ID and return the appropriate uptake bound.
fn uptake_bound_for_metabolite(met_id: &str, bounds: &MediumBounds) -> f64 {
    let lower = met_id.to_lowercase();
    let clean = lower.strip_prefix("m_").unwrap_or(&lower);
    let base = clean
        .trim_end_matches("_e").trim_end_matches("_e0")
        .trim_end_matches("_c").trim_end_matches("_p");

    // ── Amino acids ────────────────────────────────────────────────────
    // Standard 20 + selenocysteine + pyrrolysine + common derivatives
    const AMINO_ACIDS: &[&str] = &[
        "ala__l", "arg__l", "asn__l", "asp__l", "cys__l",
        "gln__l", "glu__l", "gly",    "his__l", "ile__l",
        "leu__l", "lys__l", "met__l", "phe__l", "pro__l",
        "ser__l", "thr__l", "trp__l", "tyr__l", "val__l",
        "ala__d", "arg__d", "asn__d", "asp__d", "cys__d",
        "gln__d", "glu__d", "his__d", "ile__d", "leu__d",
        "lys__d", "met__d", "phe__d", "pro__d", "ser__d",
        "thr__d", "trp__d", "tyr__d", "val__d",
        "orn",    "citr",   "hcys__l", "gaba",  "taurine",
        "cgly",   "5mthf",
    ];
    // SEED cpd IDs for amino acids (gapseq models)
    const AMINO_ACIDS_SEED: &[&str] = &[
        "cpd00023", "cpd00033", "cpd00035", "cpd00039", "cpd00041",
        "cpd00051", "cpd00053", "cpd00054", "cpd00060", "cpd00065",
        "cpd00066", "cpd00069", "cpd00084", "cpd00107", "cpd00117",
        "cpd00119", "cpd00129", "cpd00132", "cpd00135", "cpd00156",
        "cpd00161", "cpd00186", "cpd00268", "cpd00320", "cpd00322",
        "cpd00345", "cpd00549", "cpd00550", "cpd00567", "cpd00586",
        "cpd00587", "cpd00637", "cpd01017", "cpd19016",
    ];
    if AMINO_ACIDS.contains(&base) || AMINO_ACIDS_SEED.contains(&base) {
        return bounds.amino_acid;
    }

    // ── Nucleobases / nucleosides / nucleotides ────────────────────────
    const NUCLEOBASES: &[&str] = &[
        "ura", "gua", "ade", "cyt", "thy", "hyxn",
        "uri", "guo", "ado", "cyd", "thd", "ino",
        "xan", "xanthine",
        "amp", "gmp", "cmp", "ump", "imp",
        "adn", "gsn", "csn", "urd",
        "duri", "dguo", "dado", "dcyt", "dtmp", "dump",
    ];
    // SEED cpd IDs for nucleobases/nucleosides (gapseq models)
    const NUCLEOBASES_SEED: &[&str] = &[
        "cpd00018", "cpd00046", "cpd00091", "cpd00092", "cpd00114",
        "cpd00126", "cpd00128", "cpd00182", "cpd00207", "cpd00249",
        "cpd00298", "cpd00299", "cpd00307", "cpd00309", "cpd00311",
        "cpd00412", "cpd00654",
    ];
    if NUCLEOBASES.contains(&base) || NUCLEOBASES_SEED.contains(&base) {
        return bounds.nucleobase;
    }

    // ── Vitamins and cofactors ──────────────────────────────────────────
    const COFACTORS: &[&str] = &[
        "fol", "thf", "5mthf", "mlthf",
        "nmn", "nad", "nadp", "fmn", "fad", "ribflv",
        "cbl1", "cbl2", "adocbl", "mecbl",
        "btn", "pnto__r", "pnto__s", "4ppcys", "pap",
        "pyam5p", "pydx", "pydxn", "pydx5p",
        "thm", "thmpp",
        "pheme", "sheme", "hem",
        "mqn7", "mqn8", "mql8", "q8", "q8h2",
        "lipoate", "lipoamide",
        "spmd", "ptrc", "cadav",
    ];
    // SEED cpd IDs for vitamins/cofactors (gapseq models)
    const COFACTORS_SEED: &[&str] = &[
        "cpd00003", "cpd00006", "cpd00015", "cpd00016", "cpd00028",
        "cpd00045", "cpd00050", "cpd00056", "cpd00087", "cpd00104",
        "cpd00118", "cpd00125", "cpd00166", "cpd00215", "cpd00220",
        "cpd00263", "cpd00264", "cpd00305", "cpd00355", "cpd00393",
        "cpd00423", "cpd00493", "cpd00541", "cpd00557", "cpd00635",
        "cpd00644", "cpd02666", "cpd11606", "cpd15499", "cpd15500",
        "cpd15561",
    ];
    if COFACTORS.contains(&base) || COFACTORS_SEED.contains(&base) {
        return bounds.cofactor;
    }

    // ── Default: carbon source / other ─────────────────────────────────
    bounds.carbon_source
}

// ─────────────────────────────────────────────────────────────────────
// CSV medium file parser
// Format: compounds,name,maxFlux  (SEED compound IDs, one per line)
// ─────────────────────────────────────────────────────────────────────

/// Parse a CSV medium file with the format:
/// ```text
/// compounds,name,maxFlux
/// cpd00027,D-Glucose,0.5
/// cpd00013,NH3,5
/// ...
/// ```
///
/// Returns `(compounds, per_compound_bounds)` where:
/// - `compounds` is the set of SEED compound IDs (base IDs, no compartment suffix)
/// - `per_compound_bounds` maps each base SEED ID to its max uptake flux (mmol/gDW/h)
///
/// The `maxFlux` is the **import** upper bound.  Zero-flux rows are skipped.
pub fn parse_medium_csv(path: &str) -> anyhow::Result<(HashSet<String>, HashMap<String, f64>)> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("Cannot read medium CSV '{}': {}", path, e))?;
    let mut compounds: HashSet<String> = HashSet::new();
    let mut per_bounds: HashMap<String, f64> = HashMap::new();
    for (i, line) in content.lines().enumerate() {
        if i == 0 { continue; }               // skip header
        let line = line.trim();
        if line.is_empty() { continue; }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() < 3 { continue; }
        let cpd_id = fields[0].trim().to_string();
        if cpd_id.is_empty() { continue; }
        let max_flux: f64 = fields[2].trim().parse().unwrap_or(0.0);
        if max_flux > 0.0 {
            compounds.insert(cpd_id.clone());
            per_bounds.insert(cpd_id, max_flux);
        }
    }
    if compounds.is_empty() {
        anyhow::bail!("Medium CSV '{}' contains no valid compounds", path);
    }
    eprintln!("  Medium CSV: {} compounds loaded from '{}'", compounds.len(), path);
    Ok((compounds, per_bounds))
}

// ─────────────────────────────────────────────────────────────────────
// Apply medium to model
// Equivalent to COBRApy's `model.medium = {...}` property setter
// ─────────────────────────────────────────────────────────────────────

/// Apply medium constraints to the model, replicating COBRApy exactly.
///
/// ## Algorithm (matches COBRApy `Model.medium` setter)
///
/// ```python
/// # COBRApy source (cobra/core/model.py):
/// def set_active_bound(reaction, bound):
///     if reaction.reactants:          # stoich < 0
///         reaction.lower_bound = -bound
///     elif reaction.products:         # stoich > 0
///         reaction.upper_bound = bound
///
/// for rxn in self.exchanges:          # Step 1: CLOSE ALL
///     set_active_bound(rxn, 0)
/// for rxn_id, bound in medium.items():# Step 2: OPEN MEDIUM
///     set_active_bound(rxn, bound)
/// ```
///
/// ## Why this matters
///
/// Without Step 1, exchange reactions retain their SBML default bounds
/// (typically lb=-1000, ub=1000), meaning **ALL nutrients are available**.
/// This completely invalidates the defined medium.
///
/// Returns `(n_exchanges, n_opened, n_closed)` for diagnostics.
///
/// `per_compound_bounds` — when `Some`, overrides the tiered `bounds` lookup
/// for each compound that appears in the map (keyed by SEED base ID such as
/// `cpd00027`).  This is used when loading a CSV medium that specifies an
/// explicit `maxFlux` per compound.  Compounds absent from the map fall back
/// to the tiered `bounds` classification.
pub fn apply_medium(
    model: &mut MetabolicModel,
    medium_compounds: &HashSet<String>,
    bounds: &MediumBounds,
    per_compound_bounds: Option<&HashMap<String, f64>>,
) -> (usize, usize, usize) {
    let expanded = expand_medium_compounds(medium_compounds);
    let exchange_ids = find_exchange_reactions(model);
    let n_exchanges = exchange_ids.len();
    let mut n_opened = 0usize;
    let mut n_closed = 0usize;
    let mut opened_names: Vec<String> = Vec::new();

    for rid in &exchange_ids {
        if let Some(rxn) = model.reactions.get_mut(rid) {
            if rxn.metabolites.is_empty() {
                continue;
            }
            let met_id = rxn.metabolites[0].0.clone();
            let stoich = rxn.metabolites[0].1;
            let in_medium = metabolite_in_medium(&met_id, &expanded);

            // ── Direction handling ──
            //
            // Typical exchange: "met_e →" (met is reactant, stoich < 0)
            //   Forward  = export  (positive flux)
            //   Reverse  = import  (negative flux → lower_bound)
            //
            // Atypical:  "→ met_e" (met is product, stoich > 0)
            //   Forward  = import  (positive flux → upper_bound)
            //   Reverse  = export
            //
            // COBRApy only touches the IMPORT-direction bound;
            // the EXPORT-direction bound stays at its SBML value.

            // Per-compound CSV bounds take priority over tiered classification.
            // Strip M_ prefix and compartment suffix to get the base SEED ID.
            let uptake_bound = if let Some(per) = per_compound_bounds {
                let lower = met_id.to_lowercase();
                let clean = lower.strip_prefix("m_").unwrap_or(&lower);
                let base_id = clean
                    .trim_end_matches("_e0").trim_end_matches("_e")
                    .trim_end_matches("_c0").trim_end_matches("_c")
                    .trim_end_matches("_p0").trim_end_matches("_p");
                per.get(base_id).copied()
                    .unwrap_or_else(|| uptake_bound_for_metabolite(&met_id, bounds))
            } else {
                uptake_bound_for_metabolite(&met_id, bounds)
            };

            // Check if this is a cofactor/vitamin exchange that the model
            // pre-opened in SBML (lb < 0). If not in medium but pre-opened
            // as a cofactor, preserve the SBML lb instead of closing — these
            // trace micronutrients are present in every realistic growth medium
            // and their exclusion causes gapseq models to show zero growth.
            let sbml_lb = rxn.lower_bound;
            let sbml_ub = rxn.upper_bound;
            let is_cofactor = is_cofactor_metabolite(&met_id);

            if stoich < 0.0 {
                // Metabolite is reactant: import = reverse ⇒ lower_bound
                if in_medium {
                    rxn.lower_bound = -uptake_bound;
                    opened_names.push(format!("{}({}, ub={:.1})", rid, met_id, uptake_bound));
                    n_opened += 1;
                } else if is_cofactor && sbml_lb < -1e-9 {
                    // Preserve SBML-defined cofactor uptake — do not close
                    n_opened += 1; // count as "open" for diagnostics
                } else {
                    rxn.lower_bound = 0.0; // close import
                    n_closed += 1;
                }
                // upper_bound (export direction) left unchanged
            } else {
                // Metabolite is product: import = forward ⇒ upper_bound
                if in_medium {
                    rxn.upper_bound = uptake_bound;
                    opened_names.push(format!("{}({}, ub={:.1})", rid, met_id, uptake_bound));
                    n_opened += 1;
                } else if is_cofactor && sbml_ub > 1e-9 {
                    // Preserve SBML-defined cofactor uptake — do not close
                    n_opened += 1;
                } else {
                    rxn.upper_bound = 0.0; // close import
                    n_closed += 1;
                }
                // lower_bound (export direction) left unchanged
            }
            let _ = (sbml_lb, sbml_ub); // silence unused warnings
        }
    }

    eprintln!(
        "  Medium: {} exchanges found, {} opened, {} closed ({} compounds → {} expanded IDs)",
        n_exchanges,
        n_opened,
        n_closed,
        medium_compounds.len(),
        expanded.len(),
    );
    if n_opened <= 30 {
        eprintln!("  Medium opened: {}", opened_names.join(", "));
    }

    (n_exchanges, n_opened, n_closed)
}

// ─────────────────────────────────────────────────────────────────────
// Diagnostics — use to compare with COBRApy output
// ─────────────────────────────────────────────────────────────────────

/// Print per-exchange-reaction state for debugging.
///
/// Compare with COBRApy:
/// ```python
/// for r in model.exchanges:
///     met = list(r.metabolites.keys())[0]
///     print(f"{r.id}\t{met.id}\t{r.lower_bound}\t{r.upper_bound}")
/// ```
#[allow(dead_code)]
pub fn dump_exchange_diagnostics(
    model: &MetabolicModel,
    medium_compounds: &HashSet<String>,
) {
    let expanded = expand_medium_compounds(medium_compounds);
    let exchange_ids = find_exchange_reactions(model);

    eprintln!(
        "  ┌── Exchange diagnostics ({} reactions) ──",
        exchange_ids.len()
    );

    let mut n_import_open = 0usize;
    let mut n_import_closed = 0usize;

    for rid in &exchange_ids {
        if let Some(rxn) = model.reactions.get(rid) {
            if rxn.metabolites.is_empty() {
                continue;
            }
            let met_id = &rxn.metabolites[0].0;
            let stoich = rxn.metabolites[0].1;
            let in_med = metabolite_in_medium(met_id, &expanded);

            let import_open = if stoich < 0.0 {
                rxn.lower_bound < -1e-9
            } else {
                rxn.upper_bound > 1e-9
            };

            if import_open {
                n_import_open += 1;
            } else {
                n_import_closed += 1;
            }

            let flag = match (in_med, import_open) {
                (true, true) => "  OK  ",
                (false, false) => "  OK  ",
                (true, false) => "⚠SHUT ",
                (false, true) => "⚠LEAK ",
            };

            eprintln!(
                "  │ {} {:40} met={:30} s={:5.1} lb={:8.1} ub={:8.1} med={}",
                flag, rid, met_id, stoich, rxn.lower_bound, rxn.upper_bound, in_med,
            );
        }
    }

    eprintln!(
        "  └── Summary: {} import-open, {} import-closed ──",
        n_import_open, n_import_closed,
    );
}
// ───────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Gene, MetabolicModel, Metabolite, Reaction};
    use std::collections::HashSet;

    /// Build a minimal model with a single extracellular metabolite and one
    /// exchange reaction (stoich −1 on the metabolite, the typical "(met) →"
    /// form).  Lets each test exercise apply_medium on a single exchange.
    fn make_model_with_one_ex(
        met_id: &str, ex_lb: f64, ex_ub: f64,
    ) -> MetabolicModel {
        let mut mets = HashMap::new();
        mets.insert(met_id.to_string(), Metabolite {
            id: met_id.into(), name: met_id.into(),
            compartment: Some("e0".into()), boundary: false, formula: None,
        });
        let mut rxns = HashMap::new();
        let ex_id = format!("EX_{}", met_id);
        rxns.insert(ex_id.clone(), Reaction {
            id: ex_id, name: format!("exchange {}", met_id),
            metabolites: vec![(met_id.into(), -1.0)],
            reversible: true, lower_bound: ex_lb, upper_bound: ex_ub,
            gene_reaction_rule: String::new(), gene_ids: vec![],
        });
        MetabolicModel {
            id: "test".into(), name: None,
            metabolites: mets, reactions: rxns,
            biomass_reaction: None,
            genes: HashMap::<String, Gene>::new(),
        }
    }

    #[test]
    fn test_cofactor_pre_opened_uptake_is_preserved() {
        // cpd00220 = riboflavin (a known SEED-cpd cofactor).  Suppose the SBML
        // pre-opened it with lb = −0.1.  The user-declared medium does NOT
        // include it.  apply_medium MUST keep the import bound at lb = −0.1
        // rather than closing it to 0, because vitamins / cofactors are
        // present in every realistic growth medium even when not listed.
        let mut m = make_model_with_one_ex("cpd00220_e0", -0.1, 1000.0);
        let medium: HashSet<String> = HashSet::new(); // empty
        let bounds = MediumBounds::default();
        apply_medium(&mut m, &medium, &bounds, None);
        let ex = m.reactions.get("EX_cpd00220_e0").unwrap();
        assert!(
            (ex.lower_bound - (-0.1)).abs() < 1e-9,
            "cofactor uptake bound must be preserved (got lb={})", ex.lower_bound
        );
    }

    #[test]
    fn test_non_cofactor_not_in_medium_is_closed() {
        // cpd00027 = D-glucose (a carbon source, NOT a cofactor).  When NOT
        // in the medium, apply_medium MUST close uptake (lb → 0) regardless
        // of the SBML pre-open bound.
        let mut m = make_model_with_one_ex("cpd00027_e0", -10.0, 1000.0);
        let medium: HashSet<String> = HashSet::new();
        let bounds = MediumBounds::default();
        apply_medium(&mut m, &medium, &bounds, None);
        let ex = m.reactions.get("EX_cpd00027_e0").unwrap();
        assert!(
            ex.lower_bound.abs() < 1e-9,
            "non-cofactor uptake must be closed (got lb={})", ex.lower_bound
        );
    }

    #[test]
    fn test_cofactor_in_medium_uses_tier_bound() {
        // When a cofactor IS in the user-declared medium, apply_medium should
        // open it at the cofactor TIER bound (default 0.1), not at whatever
        // value the SBML had pre-opened.
        let mut m = make_model_with_one_ex("cpd00220_e0", -0.001, 1000.0);
        let mut medium = HashSet::new();
        medium.insert("cpd00220_e0".to_string());
        let bounds = MediumBounds::default(); // cofactor = 0.1
        apply_medium(&mut m, &medium, &bounds, None);
        let ex = m.reactions.get("EX_cpd00220_e0").unwrap();
        assert!(
            (ex.lower_bound - (-bounds.cofactor)).abs() < 1e-9,
            "cofactor in medium must use tier bound (got lb={}, expected {})",
            ex.lower_bound, -bounds.cofactor
        );
    }

    #[test]
    fn test_cofactor_not_pre_opened_stays_closed() {
        // Symmetric guard: if a cofactor is NOT in the medium AND was NOT
        // pre-opened in SBML (lb = 0), apply_medium must NOT spontaneously
        // open it. Otherwise downstream FBA leaks free cofactor uptake.
        let mut m = make_model_with_one_ex("cpd00220_e0", 0.0, 1000.0);
        let medium: HashSet<String> = HashSet::new();
        let bounds = MediumBounds::default();
        apply_medium(&mut m, &medium, &bounds, None);
        let ex = m.reactions.get("EX_cpd00220_e0").unwrap();
        assert!(
            ex.lower_bound.abs() < 1e-9,
            "cofactor NOT pre-opened must stay closed (got lb={})", ex.lower_bound
        );
    }
}
