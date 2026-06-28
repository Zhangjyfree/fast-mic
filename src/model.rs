use serde::{Serialize, Deserialize};
use std::collections::HashMap;

// ─────────────────────────────────────────────────────────────────────────────
// Shared compartment + exchange constants
// Single source of truth used by medium.rs, community.rs, and cobra.rs.
// ─────────────────────────────────────────────────────────────────────────────

/// Canonical list of external compartment identifiers.
/// Merges COBRApy's `compartment_shortlist["e"]` (medium.rs) and the
/// community-model recogniser (community.rs) into one list so both modules
/// agree on what "extracellular" means.
pub const KNOWN_EXTERNAL_COMPARTMENTS: &[&str] = &[
    "e",
    "e0",           // gapseq / ModelSEED
    "C_e",          // CarveMe
    "c_e",
    "ex",
    "extr",
    "extra",
    "extracellular",
    "extraorganism",
    "out",
    "extracellular space",
    "extra organism",
    "extra cellular",
    "extra-organism",
    "external",
    "external medium",
];

/// Return true when `comp` is a recognised external compartment identifier.
/// Matching is case-insensitive and tolerates surrounding brackets.
pub fn is_extracellular_compartment(comp: &str) -> bool {
    let c = comp.trim().trim_start_matches('[').trim_end_matches(']');
    KNOWN_EXTERNAL_COMPARTMENTS.iter().any(|s| s.eq_ignore_ascii_case(c))
}

/// Reaction ID substrings that disqualify a boundary reaction from being an
/// exchange reaction.  From COBRApy `annotations.py` → `excludes["exchange"]`.
/// Matching is case-sensitive, matching COBRApy behaviour.
pub const EXCHANGE_EXCLUDES: &[&str] = &[
    "demand",
    "DM_",
    "biosynthesis",
    "transcription",
    "replication",
    "SN_",
    "SK_",
    "sink",
];

/// Canonical COBRApy-style exchange reaction check.
///
/// A reaction qualifies when ALL of:
///   1. Exactly 1 metabolite (boundary / sink reaction)
///   2. No excluded substring in the reaction ID
///   3. The sole metabolite is in `external_compartment`
///      — annotation (`met.compartment`) takes priority over ID suffix
///        so models with non-standard suffixes but correct SBML annotations
///        are handled correctly.
pub fn is_canonical_exchange(
    rxn: &Reaction,
    model: &MetabolicModel,
    external_compartment: &str,
) -> bool {
    if rxn.metabolites.len() != 1 {
        return false;
    }
    if EXCHANGE_EXCLUDES.iter().any(|ex| rxn.id.contains(ex)) {
        return false;
    }
    let met_id = &rxn.metabolites[0].0;
    // Annotation-first lookup — the authoritative source when present.
    if let Some(met) = model.metabolites.get(met_id) {
        if let Some(ref comp) = met.compartment {
            return comp == external_compartment;
        }
    }
    // Fallback: check ID suffix when annotation is absent.
    let ext_suffix = format!("_{}", external_compartment);
    met_id.ends_with(&ext_suffix) || met_id.ends_with("_e")
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetabolicModel {
    pub id: String,
    pub name: Option<String>,
    #[serde(skip)]
    pub metabolites: HashMap<String, Metabolite>,
    #[serde(skip)]
    pub reactions: HashMap<String, Reaction>,
    pub biomass_reaction: Option<String>,
    #[serde(skip)]
    pub genes: HashMap<String, Gene>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Metabolite {
    pub id: String,
    pub name: String,
    pub compartment: Option<String>,
    /// If true, this species is NOT subject to steady-state mass balance (Sv=0).
    /// Corresponds to SBML `boundaryCondition="true"`.
    #[serde(default)]
    pub boundary: bool,
    /// Chemical formula from SBML FBC v2 `fbc:chemicalFormula` (e.g. `"C6H12O6"`).
    /// `None` when the model omits the attribute. Downstream code should
    /// treat `None` as "unknown" rather than failing.
    #[serde(default)]
    pub formula: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Reaction {
    pub id: String,
    pub name: String,
    pub metabolites: Vec<(String, f64)>,
    pub reversible: bool,
    pub lower_bound: f64,
    pub upper_bound: f64,
    pub gene_reaction_rule: String,
    pub gene_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: String,
    pub label: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum InteractionType {
    Mutualism,
    Commensalism,
    Parasitism,
    Competition,
    Neutral,
    Amensalism,
}

impl std::fmt::Display for InteractionType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InteractionType::Mutualism => write!(f, "mutualism"),
            InteractionType::Commensalism => write!(f, "commensalism"),
            InteractionType::Parasitism => write!(f, "parasitism"),
            InteractionType::Competition => write!(f, "competition"),
            InteractionType::Neutral => write!(f, "neutral"),
            InteractionType::Amensalism => write!(f, "amensalism"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetaboliteExchange {
    pub metabolite_id: String,
    pub flux: f64,
    pub donor_genes: Vec<String>,
    pub receiver_genes: Vec<String>,
    pub is_essential: bool,
    /// **Detection tier flag.**
    /// - `false` → direct cross-feeding: detected by sign-based mass
    ///   balance on the merged stoichiometry (donor produces, receiver
    ///   consumes the same extracellular metabolite at non-zero flux).
    ///   Strong evidence — the metabolite physically flowed between
    ///   organisms in the LP solution.
    /// - `true`  → inferred (facilitation) cross-feeding: the receiver
    ///   benefits in co-culture (Δμ > 0) but no direct sign-based
    ///   donor/receiver contribution was found. The fallback heuristic
    ///   then attributes the benefit to a candidate substrate the
    ///   donor either secretes or releases by ceasing to consume.
    ///   Weaker evidence — treat as a hypothesis, not a measurement.
    #[serde(default)]
    pub inferred: bool,
    /// **Low-confidence cross-feed flag.**
    /// Set `true` when the metabolite is one of the "sometimes-real-but-
    /// usually-artifact" species (pyruvate, lactate, ethanol, glycerol).
    /// We intentionally don't FILTER these — they CAN be genuine in
    /// some anaerobic contexts — but downstream tools should render them
    /// differently (half-alpha, dashed, parenthesised) so that reviewers
    /// see we considered the artifact risk. Orthogonal to `inferred`.
    #[serde(default)]
    pub low_confidence: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SharedResource {
    pub metabolite_id: String,
    pub uptake_a: f64,
    pub uptake_b: f64,
    pub total_available: f64,
    pub competition_ratio: f64,
    pub genes_a: Vec<String>,
    pub genes_b: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetFluxRecord {
    /// The reaction ID specified by the user (e.g., "EX_ppa_e").
    pub reaction_id: String,
    /// Flux when species A is grown alone (mmol/gDW/h).
    /// For exchange reactions: positive = secretion, negative = uptake.
    pub alone_a: f64,
    /// Flux when species B is grown alone.
    pub alone_b: f64,
    /// Net community exchange flux in co-culture
    /// (the shared EX_* reaction value — sign convention same as alone_*).
    pub co_total: f64,
    /// Species A's contribution to the extracellular metabolite in co-culture,
    /// computed as Σ(stoichiometry × flux) over all A-suffixed reactions
    /// touching the target metabolite. Positive = net production by A.
    pub co_a_contribution: f64,
    /// Species B's contribution (same convention).
    pub co_b_contribution: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseResult {
    pub species_a: String,
    pub species_b: String,
    pub growth_a_alone: f64,
    pub growth_b_alone: f64,
    pub growth_a_co: f64,
    pub growth_b_co: f64,
    pub mu: f64,
    pub reciprocity_index: f64,
    pub benefit_a: f64,
    pub benefit_b: f64,
    pub net_interaction: f64,
    pub interaction_type: InteractionType,
    pub cross_feeding_score: f64,

    /// Fraction of total cross-feeding flux that has ≥1 real gene
    /// attribution on either donor or receiver side. Range [0.0, 1.0]:
    ///   1.0 = every cf-flux mmol is GPR-supported,
    ///   0.0 = no cf-flux has any gene attribution,
    ///   intermediate values = partial coverage.
    /// `NaN` when the pair has no detected cross-feeding (i.e. the
    /// metric is undefined for this pair). Downstream may bin this
    /// continuous value as desired (e.g. ≥0.9 → "high-confidence").
    pub gene_supported_fraction: f64,

    pub a_to_b_exchanges: Vec<MetaboliteExchange>,
    pub b_to_a_exchanges: Vec<MetaboliteExchange>,
    pub shared_uptakes: Vec<SharedResource>,
    #[serde(default)]
    pub target_fluxes: Vec<TargetFluxRecord>,

    // ── Continuous niche / competition metrics (MIND-inspired) ────────────
    /// Jaccard overlap of monoculture-uptake substrates between A and B,
    /// approximated from the pairwise co-culture: |shared| / |shared ∪
    /// directional exchanges|. Range [0, 1]. Higher = more substrate
    /// niche overlap → expected stronger competition. `None` if uptake
    /// sets cannot be reliably inferred (e.g. one species non-viable).
    #[serde(default)]
    pub niche_overlap_jaccard: Option<f64>,

    /// Flux-weighted niche overlap. Sum over shared substrates of
    /// `min(uptake_a, uptake_b)` divided by the total uptake flux of
    /// the lower-uptake species (incl. cross-feeding). Range [0, 1];
    /// emphasises substrates where both species actually invest
    /// significant uptake capacity.
    #[serde(default)]
    pub niche_overlap_weighted: Option<f64>,

    /// Competition intensity (Ruzicka / weighted-Jaccard similarity)
    /// computed over the **union** of extracellular substrates utilised
    /// by either species in co-culture:
    ///
    /// ```text
    ///     Σ min(u_a, u_b) / Σ max(u_a, u_b)
    /// ```
    ///
    /// where u_x = 0 for substrates the species does not draw on.
    /// Bounded in [0, 1] and dimensionless, so directly comparable
    /// across pairs of any biomass / flux scale:
    ///   0 → niches are disjoint (no shared substrate use),
    ///   1 → identical substrate use at identical rates (maximal
    ///       competition).
    /// Sensitive to **both** breadth (disjoint substrates inflate the
    /// denominator but not the numerator) **and** rate mismatch on
    /// shared substrates. This is the flux-weighted analogue of
    /// `niche_overlap_jaccard` and the symmetric counterpart of the
    /// asymmetric `niche_overlap_weighted`.
    #[serde(default)]
    pub competition_intensity: f64,
}