# fast-mic

**Fast Metabolic Interaction Calculator** — A high-performance Rust tool for pairwise pFBA-based microbial community analysis. Compute cross-feeding, competition, and interaction-type predictions from genome-scale metabolic models (GEMs) at scale.

> **Associated study**: fast-mic was developed for the study *"Fast-mic: a scalable tool for exhaustive pairwise interaction typing of genome-scale metabolic models reveals that carbon quality shapes probiotic–microbiome cooperation across a prebiotic gradient"* — an exhaustive pairwise interaction screen of six *Akkermansia* strains (3 species, mucin specialist) and ten *Lactobacillus*-group strains (9 species, metabolic generalist) against the gut (UHGG) community across a 10-level prebiotic gradient (L0–L9).

---

## Installation

```bash
git clone https://github.com/Zhangjyfree/fast-mic.git
cd fast-mic
cargo build --release
# Binaries: ./target/release/{fast-mic, bench-single-fba, bench-cff-deviation}
```

Requires Rust ≥ 1.75 and a C compiler for the bundled HiGHS solver.

---

## Quick start

```bash
# All-vs-all pairwise within one model set
fast-mic \
  --medium-name WesternDiet \
  --media-db media/media_db.tsv \
  --compounds-tsv media/compounds.tsv \
  -o pairwise.tsv \
  models/*.xml

# Cross-group: each Akkermansia × each gut commensal
fast-mic \
  --group1 akk_strains/ \
  --group2 commensals/ \
  --medium-file media/western_diet_mucin_gapseq.csv \
  -o akk_vs_gut.tsv --full-tsv akk_vs_gut_full.tsv \
  --target-reaction EX_ppa_e,EX_ac_e,EX_but_e \
  --threads 0
```

---

## Command-line reference

### Input

| Flag | Description |
|---|---|
| `<files>...` | Positional SBML model paths; forms all-vs-all pairs. |
| `--group1 <DIR>` | Directory of `.xml`/`.sbml` files (group 1). |
| `--group2 <DIR>` | Directory of `.xml`/`.sbml` files (group 2). Computes cross-group pairs with `--group1`. |
| `--medium-name <NAME>` | Named medium from `--media-db`. Mutually exclusive with `--medium-file`. |
| `--medium-file <FILE>` | CSV medium file with SEED compound IDs and per-compound `maxFlux`. For gapseq models. |
| `--media-db <FILE>` | Medium definition TSV (`medium`, `description`, `compound`, `name`). Default: `media_db.tsv`. |
| `--compounds-tsv <FILE>` | ModelSEED `compounds.tsv` for BiGG→SEED translation. Default: `compounds.tsv`. |
| `--pair-filter <FILE>` | 2-column TSV restricting which `(species_a, species_b)` pairs to compute. |

### Output

| Flag | Description |
|---|---|
| `-o, --output <FILE>` | Compact pairwise TSV. Default: `output.tsv`. |
| `--full-tsv <FILE>` | Verbose TSV with per-metabolite cross-feeding details and gene attributions. |
| `--json <FILE>` | JSON dump of all pairwise results. |
| `-v, --verbose` | Per-pair details to stderr. |
| `--summary` | Suppress per-pair output; print final summary only. |

### Medium uptake limits

| Flag | Default | Description |
|---|---|---|
| `--medium-uptake-limit` | 10.0 | Max uptake rate (mmol/gDW/h) for carbon-source compounds. |

Tiered limits for amino acids (1.0), nucleobases/nucleosides (0.5), and cofactors (0.1) are applied automatically based on compound classification.

### Target-reaction tracking

| Flag | Description |
|---|---|
| `--target-reaction R1,R2,...` | Track flux of one or more reactions (e.g. `EX_ac_e,EX_ppa_e,EX_but_e` for SCFAs). Outputs five columns per reaction: `alone_a`, `alone_b`, `co_total`, `co_a`, `co_b`. |

### LP tolerances

| Flag | Default | Description |
|---|---|---|
| `--lock-tol` | 1e-5 | Tolerance for pinning fluxes in CFF/pFBA lock constraints (`v ∈ [v* ± tol]`). 100 × HiGHS feasibility tolerance. Drop to 1e-7 only for exact single-species reproduction. |

### Performance

| Flag | Default | Description |
|---|---|---|
| `--threads N` | 0 | Worker threads. 0 = all cores; 1 = serial. |
| `--cache-monoculture` | true | Pre-compute monoculture pFBA per model and reuse across pairs (~30-50 % speedup). |

### Co-culture objective

| Flag | Default | Description |
|---|---|---|
| `--fixed-ratio` | off | Use a fixed-ratio co-culture objective — total community biomass maximised subject to μ_A/μ_B pinned to the monoculture ratio (μ_A^co/μ_B^co = μ_A^alone/μ_B^alone) — instead of the default lexicographic max-min allocation. Intended for objective-sensitivity analysis. |

---

## Output schemas

### Compact TSV (`-o`)

Each row is one species pair.

| Column | Description |
|---|---|
| `species_a`, `species_b` | Model IDs. |
| `growth_a_alone`, `growth_b_alone` | Monoculture growth rates (h⁻¹). |
| `growth_a_co`, `growth_b_co` | Co-culture growth rates. |
| `benefit_a`, `benefit_b` | `(growth_co − growth_alone) / growth_alone`. |
| `interaction_type` | `mutualism`, `commensalism`, `parasitism`, `competition`, `amensalism`, or `neutral`. |
| `gene_supported_fraction` | Fraction of cross-feeding flux attributable to annotated genes. |
| `n_exchanged_metabolites` | Number of metabolites exchanged in either direction. |
| `competition_intensity` | Σ min(uptake_a, uptake_b) over shared resources. |
| `{rxn}__alone_a/b`, `{rxn}__co_total/a/b` | Per-reaction flux columns for each `--target-reaction`. |

### Full TSV (`--full-tsv`)

Adds per-metabolite cross-feeding columns: `a_to_b_metabolites`, `a_to_b_fluxes`, `a_to_b_donor_genes`, `a_to_b_receiver_genes`, mirror columns for B→A, `a_to_b_inferred` (hypothesis-grade entries), `a_to_b_low_confidence`, and competed-resource columns. Lists are `;`-separated.

---

## Algorithm

### Single-species: FBA → CycleFreeFlux + pFBA

After standard FBA (max biomass), a single LP achieves cycle removal **and** parsimony simultaneously:

```
min  Σ |v_i|        over non-exchange, non-biomass reactions
s.t. S v = 0
     v_exch_j ∈ [v*_exch_j ± ε]              (preserve FBA exchange profile)
     v*_biomass − ε ≤ v_biomass ≤ v*_biomass + ε   (hold growth at the FBA optimum)
     lb_i ≤ v_i ≤ ub_i
```

Fixing exchange fluxes eliminates Type-III internal cycles (they carry net-zero exchange flux). Minimising Σ|v| drives every closed loop to zero. The biomass flux is **constrained to its Stage-1 FBA optimum within the lock tolerance ε** (band form, `ε = LOCK_TOL`), so the parsimony objective cannot trade growth rate for lower total flux — loop removal therefore preserves the optimal growth rate rather than only bounding it above. This is validated empirically: across 7,463 growing model–medium evaluations (1,000 UHGG models × the L0–L9 gradient) the post-CycleFreeFlux biomass deviates from the FBA optimum by at most **1.0 × 10⁻⁵ h⁻¹ (= ε)**, biologically negligible (≤ 5.2 × 10⁻⁴ relative); see [Benchmarking](#benchmarking).

**Reference**: Desouki *et al.* (2015), *CycleFreeFlux*, BMC Bioinformatics **16**:283.

### Pairwise co-culture: lexicographic max-min

Two species' models are merged into a synthetic joint model sharing a common extracellular pool, with species-specific exchanges replaced by community-level `EX_*` reactions. Joint growth is optimised under a two-phase LP:

**LP1 — Rawlsian fairness**:
```
max z   s.t.  g_A ≥ z · g_A^alone
              g_B ≥ z · g_B^alone,  S v = 0, bounds
```

**LP2 — Utilitarian productivity within the fair set**:
```
max g_A + g_B   s.t.  g_A ≥ (z* − τ) · g_A^alone
                      g_B ≥ (z* − τ) · g_B^alone,  S v = 0, bounds
```

A co-culture CycleFreeFlux pass then removes inter-species cycles.

**Reference**: Bertsimas *et al.* (2011), *The Price of Fairness*, Operations Research **59**(1):17-31.

**Alternative objective (`--fixed-ratio`)**: For objective-sensitivity analysis, the two-phase allocation can be replaced by a single LP that maximises total community biomass subject to a fixed growth ratio, `g_A / g_B = g_A^alone / g_B^alone`. On a representative subset (*L. gasseri* × UHGG at L5 / L6) the two objectives give identical mutualism fractions (33.4 % / 22.4 %) and ≥ 99 % per-pair classification agreement (all disagreements confined to the commensalism↔neutral boundary), confirming the core conclusions are robust to the choice of growth-allocation rule.

### Design principles

Three scalar constants control the entire FBA pipeline; one is user-configurable.

| Constant | Default | Role |
|---|---|---|
| `MIN_VIABLE_GROWTH` | 1×10⁻⁴ h⁻¹ | Biological viability floor (~4 doublings/day). |
| `NUMERICAL_TOL` | 1×10⁻⁶ | LP comparison tolerance (10 × HiGHS default primal tolerance). |
| `LOCK_TOL` (`--lock-tol`) | 1×10⁻⁵ | CFF/pFBA lock-constraint tolerance — the ε band on both exchange-flux locks and the biomass constraint (`v ∈ [v* ± ε]`). Configurable. |

No empirical synergy caps, flux-ratio thresholds, metabolite blacklists, or Pareto-front scans.

---

## Supported SBML features

The parser is a hand-written two-pass `quick-xml` reader. It handles the common GEM dialects (BiGG, AGORA, CarveMe, gapseq) but does **not** implement the full SBML L3 spec.

**Supported:**
- SBML L3 core: `<model>`, `<listOfSpecies>`, `<listOfReactions>`, `<listOfCompartments>`, `<listOfParameters>`
- Species attributes: `id`, `name`, `compartment`, `boundaryCondition`, `fbc:chemicalFormula`
- Reaction attributes: `id`, `name`, `reversible`, `fbc:lowerFluxBound`, `fbc:upperFluxBound`
- FBC v2: `<fbc:listOfObjectives>`, `<fbc:objective>`, `<fbc:fluxObjective>`, `<fbc:listOfGeneProducts>`, `<fbc:geneProductAssociation>` with `<fbc:and>`/`<fbc:or>` trees
- Namespace-prefixed attribute names via local-name matching

**NOT supported** (silently ignored or fallback to defaults):
- `<initialAssignment>` — bounds set via initial assignment not picked up
- `<listOfRules>` — assignment / rate rules ignored
- `<listOfEvents>` — kinetic events ignored
- FBC v1 (deprecated upstream)

If your model uses unsupported features, convert it to plain FBC v2 first (e.g. `cobrapy.io.write_sbml_model`).

---

## Media database format

`media_db.tsv` — tab-separated TSV with four columns: `medium`, `description`, `compound`, `name`.

```
medium        description          compound    name
WesternDiet   AGORA Western Diet   glc__D      D-Glucose
WesternDiet   AGORA Western Diet   ala__L      L-Alanine
WesternDiet   AGORA Western Diet   cpd00027    D-Glucose (SEED)
LB            Lysogeny broth       ala__L      L-Alanine
```

BiGG (`glc__D`) and ModelSEED (`cpd00027`) IDs are both accepted. With `--medium-name`, fast-mic translates BiGG → SEED via `--compounds-tsv` (ModelSEED `compounds.tsv`) so AGORA-style and gapseq-style models are both matched.

### CSV medium file

For gapseq models, use `--medium-file` with CSV format:

```
compounds,name,maxFlux
cpd00027,D-Glucose,10.0
cpd00035,L-Alanine,1.0
cpd00009,Phosphate,1000.0
```

### Tiered uptake limits

Compound class is auto-detected from compound ID.

| Class | Default rate | Examples |
|---|---|---|
| Carbon sources | `--medium-uptake-limit` (10.0 mmol/gDW/h) | Sugars, organic acids |
| Amino acids | 1.0 mmol/gDW/h | Standard 20 + D-forms + ornithine |
| Nucleobases / nucleosides | 0.5 mmol/gDW/h | Adenine, uracil, adenosine |
| Cofactors / vitamins | 0.1 mmol/gDW/h | Folate, B12, riboflavin |
| Inorganic ions | Unlimited (−1000) | Na⁺, K⁺, Fe²⁺, PO₄³⁻ |

---

## Examples

### Pairwise + SCFA tracking

```bash
fast-mic \
  --group1 producers/ --group2 consumers/ \
  --medium-name WesternDiet \
  --media-db media/media_db.tsv \
  --compounds-tsv media/compounds.tsv \
  --target-reaction EX_ac_e,EX_ppa_e,EX_but_e \
  -o scfa.tsv --full-tsv scfa_full.tsv \
  --threads 0
```

---

## Prebiotic-gradient media

`media/` ships the 10 gapseq medium CSVs (`gradient_L0_base_gapseq.csv` … `gradient_L9_mos_gapseq.csv`) plus `gradient_media_list.txt` (one CSV path per line), spanning a cumulative prebiotic gradient. Each level adds one prebiotic's hydrolysis products on top of the mucin-containing base (design: *Akkermansia* viable at all levels).

| Level | Prebiotic | New compounds |
|---|---|---|
| L0 | Base (mucin) | — |
| L1 | Inulin | D-Fructose, Sucrose |
| L2 | FOS | Inulobiose |
| L3 | GOS | Lactulose, Lactose, D-Galactose |
| L4 | XOS / Arabinoxylan | D-Xylose, L-Arabinose |
| L5 | Pectin | Galacturonate, L-Rhamnose |
| L6 | Resistant starch | D-Glucose, Maltose, Maltodextrin |
| L7 | β-glucan | Cellobiose |
| L8 | HMO | Lacto-N-biose |
| L9 | MOS | D-Mannose, Mannobiose |

**Run the pairwise gradient screen, one level at a time**:

```bash
for L in media/gradient_L*_gapseq.csv; do
  fast-mic --group1 akk_strains/ --group2 commensals/ \
    --medium-file "$L" --threads 0 \
    -o "gradient_$(basename "$L" .csv).tsv"
done
```

For single-species throughput across **all** levels in one pass, the bench tools below accept `--media-list media/gradient_media_list.txt`.

> **Figures, supplementary tables, and the full COBRApy-comparison / reproduction pipeline** (R plotting scripts, `run_thread_scaling.sh`, the `cobrapy` cross-check, etc.) are **not** part of this crate — they live in the companion analysis repository. This repo ships only the Rust tool (`src/`) and the gradient media (`media/`).

---

## Benchmarking

A standalone binary `bench-single-fba` measures single-species FBA throughput across many models — under one medium or several at once — and supports thread-scaling benchmarks.

**Single medium** — a named medium from a TSV database, or one gapseq/SEED-format CSV:

```bash
# Named medium from a TSV database
bench-single-fba media/media_db.tsv WesternDiet \
  --model-list models.txt --threads 0 > bench_results.tsv

# A single gapseq/SEED-format CSV medium
bench-single-fba --medium-file media/gradient_L0_base_gapseq.csv \
  --model-list models.txt --threads 0 > bench_results.tsv
```

**Several CSV media at once** — pass `--media-list FILE`, a plain-text file with one medium-CSV path per line (absolute or relative). Each model is loaded **once** and evaluated under **every** medium, yielding one row per (model, medium) pair. This is the workload behind the 10-level prebiotic gradient.

```bash
# media/gradient_media_list.txt lists the 10 gradient CSVs (L0–L9), one path per line
bench-single-fba --media-list media/gradient_media_list.txt \
  --model-list models.txt --threads 0 > bench_gradient.tsv
```

Output columns: `model_id`, `n_metabolites`, `n_reactions`, `n_genes`, `biomass_rxn`, `growth_rate`, `load_time_s`, `fba_time_s`. With `--media-list`, `model_id` is suffixed with the medium label (the CSV file stem), e.g. `L_acidophilus_NCFM__gradient_L0_base_gapseq`, so each row stays uniquely keyed per (model, medium).

### Reference validation against COBRApy

fast-mic's single-species growth rates have been validated against COBRApy (HiGHS backend) on **9,950 genome–medium pairs** (1,000 UHGG models across the L0–L9 gradient; 5 COBRApy timeouts excluded): **Pearson r = 1.000, MAE = 3.12 × 10⁻⁷**, with 100 % of growing models agreeing to within 1 % relative error.

The end-to-end comparison driver (which additionally requires Python + `cobrapy` + `highspy`) and its threshold-asserting integration test are part of the companion analysis repository, not this crate.

### Loop-removal validation (post-CFF biomass deviation)

The binary `bench-cff-deviation` checks that the CycleFreeFlux + pFBA step does not reduce the growth rate below the Stage-1 FBA optimum. For every (model, medium) it reports the deviation `fba_optimal − post_cff_biomass`; with the band-form biomass constraint this is bounded by ε.

```bash
bench-cff-deviation \
  --media-list media/gradient_media_list.txt \
  --model-list models.txt --threads 0 > cff_deviation.tsv
```

Per-row TSV columns: `model_id`, `medium`, `fba_optimal`, `post_cff_biomass`, `deviation`, `viable`. The summary (stderr) reports the **max** and **mean** |deviation| over the growing (viable) evaluations and the worst case.

**Current status:** across 7,463 growing evaluations (1,000 models × L0–L9), max |deviation| = **1.0 × 10⁻⁵ h⁻¹ (= ε = `LOCK_TOL`)**, mean 9.7 × 10⁻⁶ h⁻¹ — within the LP tolerance and biologically negligible; growth rates are unchanged.

---

## Testing

```bash
# Unit tests (fast)
cargo test

# Lint and format
cargo clippy -- -D warnings
cargo fmt -- --check
```

Test coverage (unit tests embedded in `src/`):

| Module | Tests | Covers |
|---|---|---|
| `cobra` | 22 | Exchange detection, compartment inference, biomass finder, merged-model construction, interaction classification, NaN-safe sorts |
| `medium` | 4 | Cofactor pre-opened uptake preservation; non-cofactor closed; cofactor in medium uses tier bound; cofactor not pre-opened stays closed |

---

## Library API

`fast-mic` is also usable as a Rust library.

```rust
use fast_mic::{cobra, medium, sbml};

let model = sbml::parse_sbml("model.xml")?;
let medium_set = medium::expand_medium_compounds(&base_compounds);
let params = cobra::AnalysisParams::default();
let result = cobra::run_fba(&model, &medium_set, &params)?;
println!("growth = {}", result.objective_value);   // post-CycleFreeFlux biomass
println!("FBA optimum = {}", result.fba_optimal);  // Stage-1 optimum; difference = loop-removal deviation
```

| Module | Purpose |
|---|---|
| `model` | Core types: `MetabolicModel`, `Reaction`, `Metabolite` (with optional `formula`), `PairwiseResult`, `InteractionType`. Also exports the unified `KNOWN_EXTERNAL_COMPARTMENTS`, `EXCHANGE_EXCLUDES`, `is_canonical_exchange`, `is_extracellular_compartment`. |
| `sbml` | SBML parsing (FBC v2). Two-pass, single read. Extracts `fbc:chemicalFormula`. |
| `medium` | Compound matching, exchange-reaction detection (COBRApy-style), tiered uptake bounds, cofactor pre-opened preservation, BiGG → SEED translation. |
| `cobra` | FBA, CycleFreeFlux + pFBA (`run_fba`, `run_fba_locked`), pairwise co-culture (lexicographic max-min, or fixed-ratio via `AnalysisParams::fixed_ratio`), cross-feeding analysis. `FBAResult` carries both `objective_value` (post-CFF biomass) and `fba_optimal` (Stage-1 optimum); their difference is the loop-removal deviation. `AnalysisParams::lock_tol` and `::fixed_ratio` configurable. |

---

## Citing

If you use fast-mic in published work, please cite:

- **Desouki *et al.* (2015)**, *CycleFreeFlux: efficient removal of thermodynamically infeasible loops from flux distributions*, BMC Bioinformatics **16**:283 — *cycle-removal algorithm.*
- **Bertsimas, Farias & Trichakis (2011)**, *The Price of Fairness*, Operations Research **59**(1):17-31 — *lexicographic max-min co-culture formulation.*

---

## License

MIT

## Issues & contributions

Bug reports and PRs welcome at the project repository.
