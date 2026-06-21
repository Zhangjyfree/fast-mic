# fast-mic

**Fast Metabolic Interaction Calculator** — A high-performance Rust tool for pairwise pFBA-based microbial community analysis. Compute cross-feeding, competition, and interaction-type predictions from genome-scale metabolic models (GEMs) at scale.

快速代谢互作计算器 —— 基于 Rust 实现的高性能两两 pFBA 微生物互作分析工具。从基因组尺度代谢模型（GEM）批量计算交叉营养、竞争与互作类型。

> **Associated study / 关联研究**: fast-mic was developed for the study *"Genome-scale modelling indicates that carbon quality governs how generalist and specialist probiotics cooperate with the resident microbiome across a prebiotic gradient"* — an exhaustive pairwise interaction screen of six *Akkermansia* strains (3 species, mucin specialist) and ten *Lactobacillus*-group strains (9 species, metabolic generalist) against the gut (UHGG) community across a 10-level prebiotic gradient (L0–L9). / fast-mic 用于研究《基因组尺度建模揭示碳源质量如何决定泛化型与专化型益生菌在益生元梯度下与常驻菌群的合作》：在 10 级益生元梯度下，对 6 株 *Akkermansia*（3 个种，黏蛋白专化型）与 10 株 *Lactobacillus* 类益生菌（9 个种，代谢泛化型）与肠道（UHGG）菌群进行穷举式两两互作筛选。

---

## Highlights / 功能亮点

- **Pure Rust + HiGHS LP solver** — Single static binary, no Python / COBRApy runtime dependency. / 纯 Rust + HiGHS LP 求解器，单一静态二进制文件。
- **Hyperparameter-free FBA pipeline** — FBA → CycleFreeFlux (Desouki *et al.* 2015) → parsimonious FBA in one LP. No empirical synergy caps, no flux-ratio thresholds. / 无超参数 FBA 流程，三步合一 LP，无经验性协同上限或通量比阈值。
- **Lexicographic max-min pairwise co-culture** — Rawlsian fairness (max-min growth ratio) followed by utilitarian total-biomass maximisation. Two LPs, unique optimum, no Pareto-front scan. / 两两共培养字典序最大最小优化：先公平，后效率，唯一最优。
- **Massively parallel** — Rayon-driven per-pair parallelism, near-linear scaling to all cores. Monoculture pFBA cache delivers ~30-50 % speedup in cross-group mode. / Rayon 并行，近线性扩展，单培养 pFBA 缓存再加速 30-50 %。
- **Dual medium support** — Named medium from TSV database (BiGG IDs, auto-translated to SEED via `--compounds-tsv`) or explicit CSV file (SEED IDs + per-compound `maxFlux`, for gapseq models). / 双重培养基支持。
- **Validated against COBRApy** — `cargo test -- --ignored` enforces Pearson r ≥ 0.999, MAE ≤ 1e-3. Current cross-tool agreement on 9,950 genome–medium pairs (1,000 UHGG models × the L0–L9 gradient): r = 1.000, MAE = 3.12 × 10⁻⁷. / 与 COBRApy 在 9,950 个基因组×培养基组合上对照验证：r = 1.000、MAE = 3.12 × 10⁻⁷。
- **~221 × faster than COBRApy** on the same LP problem, single-threaded (~49 vs 0.22 genome–medium evaluations s⁻¹); peaking around ~333 evaluations s⁻¹ at 12 threads on the 1,000-model corpus. / 单线程比 COBRApy 快约 221 倍，12 线程峰值约 333 evaluations/s。
- **Cross-group analysis** — Group-1 × group-2 mode (e.g. probiotic panel × resident community) with candidate partner ranking and metabolite-flow Sankey. / 跨组分析：第一组 × 第二组（如益生菌 panel × 常驻菌群）配对，候选合作菌排名与代谢流 Sankey 图。
- **Prebiotic gradient pipeline** — 10-level cumulative medium gradient (L0–L9) covering inulin, FOS, GOS, XOS, pectin, resistant starch, β-glucan, HMO, MOS. / 10 级益生元梯度培养基流水线。
- **Interaction typing & metrics** — Per-pair relative benefit (β), competition intensity, cross-feeding score, gene-supported fraction, and six-class interaction typing (mutualism / competition / commensalism / parasitism / amensalism / neutral). / 逐对相对收益（β）、竞争强度、交叉营养得分、基因支持比例与六类互作分型。

---

## Installation / 安装

```bash
git clone https://github.com/your-org/fast-mic.git
cd fast-mic
cargo build --release
# Binary: ./target/release/fast-mic
```

Requires Rust ≥ 1.75 and a C compiler for the bundled HiGHS solver. / 需要 Rust ≥ 1.75 与 C 编译器（用于内置 HiGHS）。

---

## Quick start / 快速开始

```bash
# All-vs-all pairwise within one model set / 单组模型两两互作
fast-mic \
  --medium-name WesternDiet \
  --media-db media/media_db.tsv \
  --compounds-tsv media/compounds.tsv \
  -o pairwise.tsv \
  models/*.xml

# Cross-group: each Akkermansia × each gut commensal
# 跨组：每个 Akkermansia 模型 × 每个肠道共生菌
fast-mic \
  --group1 akk_strains/ \
  --group2 commensals/ \
  --medium-file media/western_diet_mucin_gapseq.csv \
  -o akk_vs_gut.tsv --full-tsv akk_vs_gut_full.tsv \
  --target-reaction EX_ppa_e,EX_ac_e,EX_but_e \
  --threads 0
```

---

## Command-line reference / 命令行参考

### Input / 输入

| Flag | Description / 说明 |
|---|---|
| `<files>...` | Positional SBML model paths; forms all-vs-all pairs. / 位置参数 SBML 模型路径，自动两两配对。 |
| `--group1 <DIR>` | Directory of `.xml`/`.sbml` files (group 1). / 第一组模型目录。 |
| `--group2 <DIR>` | Directory of `.xml`/`.sbml` files (group 2). Computes cross-group pairs with `--group1`. / 第二组目录，与 `--group1` 配合跨组配对。 |
| `--medium-name <NAME>` | Named medium from `--media-db`. Mutually exclusive with `--medium-file`. / 命名培养基。 |
| `--medium-file <FILE>` | CSV medium file with SEED compound IDs and per-compound `maxFlux`. For gapseq models. / SEED 格式 CSV 培养基，适用于 gapseq 模型。 |
| `--media-db <FILE>` | Medium definition TSV (`medium`, `description`, `compound`, `name`). Default: `media_db.tsv`. / 培养基定义 TSV。 |
| `--compounds-tsv <FILE>` | ModelSEED `compounds.tsv` for BiGG→SEED translation. Default: `compounds.tsv`. / 用于 BiGG→SEED ID 转换。 |
| `--pair-filter <FILE>` | 2-column TSV restricting which `(species_a, species_b)` pairs to compute. / 限制计算范围的 TSV。 |

### Output / 输出

| Flag | Description / 说明 |
|---|---|
| `-o, --output <FILE>` | Compact pairwise TSV. Default: `output.tsv`. / 紧凑两两互作 TSV。 |
| `--full-tsv <FILE>` | Verbose TSV with per-metabolite cross-feeding details and gene attributions. / 详细 TSV，含逐代谢物交叉营养与基因归因。 |
| `--json <FILE>` | JSON dump of all pairwise results. / 所有配对结果的 JSON 输出。 |
| `-v, --verbose` | Per-pair details to stderr. / 在 stderr 输出逐对详情。 |
| `--summary` | Suppress per-pair output; print final summary only. / 仅输出最终汇总。 |

### Medium uptake limits / 培养基摄取限制

| Flag | Default | Description / 说明 |
|---|---|---|
| `--medium-uptake-limit` | 10.0 | Max uptake rate (mmol/gDW/h) for carbon-source compounds. / 碳源最大摄取速率。 |

Tiered limits for amino acids (1.0), nucleobases/nucleosides (0.5), and cofactors (0.1) are applied automatically based on compound classification. / 氨基酸、核碱基/核苷、辅因子的分级限制自动应用。

### Target-reaction tracking / 目标反应追踪

| Flag | Description / 说明 |
|---|---|
| `--target-reaction R1,R2,...` | Track flux of one or more reactions (e.g. `EX_ac_e,EX_ppa_e,EX_but_e` for SCFAs). Outputs five columns per reaction: `alone_a`, `alone_b`, `co_total`, `co_a`, `co_b`. / 追踪反应的通量，每个反应输出五列。 |

### LP tolerances / LP 数值公差

| Flag | Default | Description / 说明 |
|---|---|---|
| `--lock-tol` | 1e-5 | Tolerance for pinning fluxes in CFF/pFBA lock constraints (`v ∈ [v* ± tol]`). 100 × HiGHS feasibility tolerance. Drop to 1e-7 only for exact single-species reproduction. / CFF/pFBA 锁约束公差。 |

### Performance / 性能

| Flag | Default | Description / 说明 |
|---|---|---|
| `--threads N` | 0 | Worker threads. 0 = all cores; 1 = serial. / 工作线程数。 |
| `--cache-monoculture` | true | Pre-compute monoculture pFBA per model and reuse across pairs (~30-50 % speedup). / 单培养 pFBA 缓存。 |

### Co-culture objective / 共培养目标函数

| Flag | Default | Description / 说明 |
|---|---|---|
| `--fixed-ratio` | off | Use a fixed-ratio co-culture objective — total community biomass maximised subject to μ_A/μ_B pinned to the monoculture ratio (μ_A^co/μ_B^co = μ_A^alone/μ_B^alone) — instead of the default lexicographic max-min allocation. Intended for objective-sensitivity analysis. / 使用固定比例共培养目标（共培养生长比 = 单培养比，最大化群落总生物量），替代默认字典序最大最小分配；用于目标函数敏感性分析。 |

---

## Output schemas / 输出格式

### Compact TSV (`-o`)

Each row is one species pair. / 每行为一个物种配对。

| Column | Description / 说明 |
|---|---|
| `species_a`, `species_b` | Model IDs. / 模型 ID。 |
| `growth_a_alone`, `growth_b_alone` | Monoculture growth rates (h⁻¹). / 单培养生长速率。 |
| `growth_a_co`, `growth_b_co` | Co-culture growth rates. / 共培养生长速率。 |
| `benefit_a`, `benefit_b` | `(growth_co − growth_alone) / growth_alone`. / 相对收益。 |
| `interaction_type` | `mutualism`, `commensalism`, `parasitism`, `competition`, `amensalism`, or `neutral`. / 互作类型。 |
| `gene_supported_fraction` | Fraction of cross-feeding flux attributable to annotated genes. / 有基因注释支持的比例。 |
| `n_exchanged_metabolites` | Number of metabolites exchanged in either direction. / 交换代谢物数量。 |
| `competition_intensity` | Σ min(uptake_a, uptake_b) over shared resources. / 共享资源竞争强度。 |
| `{rxn}__alone_a/b`, `{rxn}__co_total/a/b` | Per-reaction flux columns for each `--target-reaction`. / 每个目标反应的通量列。 |

### Full TSV (`--full-tsv`)

Adds per-metabolite cross-feeding columns: `a_to_b_metabolites`, `a_to_b_fluxes`, `a_to_b_donor_genes`, `a_to_b_receiver_genes`, mirror columns for B→A, `a_to_b_inferred` (hypothesis-grade entries), `a_to_b_low_confidence`, and competed-resource columns. Lists are `;`-separated.

在紧凑 TSV 基础上增加逐代谢物交叉营养列：代谢物 ID、通量、供体/受体基因、推断标记、低置信度标记、以及竞争资源列。列表以 `;` 分隔。

---

## Algorithm / 算法

### Single-species: FBA → CycleFreeFlux + pFBA

After standard FBA (max biomass), a single LP achieves cycle removal **and** parsimony simultaneously:

```
min  Σ |v_i|        over non-exchange, non-biomass reactions
s.t. S v = 0
     v_exch_j ∈ [v*_exch_j ± ε]    (preserve FBA exchange profile)
     0 ≤ v_biomass ≤ v*_biomass
     lb_i ≤ v_i ≤ ub_i
```

Fixing exchange fluxes eliminates Type-III internal cycles (they carry net-zero exchange flux). Minimising Σ|v| drives every closed loop to zero.

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

**Alternative objective (`--fixed-ratio`)**: For objective-sensitivity analysis, the two-phase allocation can be replaced by a single LP that maximises total community biomass subject to a fixed growth ratio, `g_A / g_B = g_A^alone / g_B^alone`. On a representative subset (*L. gasseri* × UHGG at L5 / L6) the two objectives give identical mutualism fractions (33.4 % / 22.4 %) and ≥ 99 % per-pair classification agreement (all disagreements confined to the commensalism↔neutral boundary), confirming the core conclusions are robust to the choice of growth-allocation rule. / 目标函数敏感性分析：可用 `--fixed-ratio` 将两阶段分配替换为"固定生长比 + 最大化总生物量"的单 LP，在代表子集上两种目标给出完全一致的互利比例与 ≥99% 的逐对分类一致率，证明核心结论不依赖分配规则。

### Design principles / 设计原则

Three scalar constants control the entire FBA pipeline; one is user-configurable.

| Constant | Default | Role / 作用 |
|---|---|---|
| `MIN_VIABLE_GROWTH` | 1×10⁻⁴ h⁻¹ | Biological viability floor (~4 doublings/day). / 生物可行性下限。 |
| `NUMERICAL_TOL` | 1×10⁻⁶ | LP comparison tolerance (10 × HiGHS default primal tolerance). / LP 比较公差。 |
| `LOCK_TOL` (`--lock-tol`) | 1×10⁻⁵ | CFF/pFBA lock-constraint tolerance. Configurable. / CFF/pFBA 锁约束公差，可配置。 |

No empirical synergy caps, flux-ratio thresholds, metabolite blacklists, or Pareto-front scans. / 无经验性协同上限、通量比阈值、代谢物黑名单或 Pareto 前沿扫描。

---

## Supported SBML features / 支持的 SBML 特性

The parser is a hand-written two-pass `quick-xml` reader. It handles the common GEM dialects (BiGG, AGORA, CarveMe, gapseq) but does **not** implement the full SBML L3 spec.

**Supported / 支持:**
- SBML L3 core: `<model>`, `<listOfSpecies>`, `<listOfReactions>`, `<listOfCompartments>`, `<listOfParameters>`
- Species attributes: `id`, `name`, `compartment`, `boundaryCondition`, `fbc:chemicalFormula`
- Reaction attributes: `id`, `name`, `reversible`, `fbc:lowerFluxBound`, `fbc:upperFluxBound`
- FBC v2: `<fbc:listOfObjectives>`, `<fbc:objective>`, `<fbc:fluxObjective>`, `<fbc:listOfGeneProducts>`, `<fbc:geneProductAssociation>` with `<fbc:and>`/`<fbc:or>` trees
- Namespace-prefixed attribute names via local-name matching

**NOT supported / 不支持** (silently ignored or fallback to defaults):
- `<initialAssignment>` — bounds set via initial assignment not picked up
- `<listOfRules>` — assignment / rate rules ignored
- `<listOfEvents>` — kinetic events ignored
- FBC v1 (deprecated upstream)

If your model uses unsupported features, convert it to plain FBC v2 first (e.g. `cobrapy.io.write_sbml_model`).

---

## Media database format / 培养基数据库格式

`media_db.tsv` — tab-separated TSV with four columns: `medium`, `description`, `compound`, `name`.

```
medium        description          compound    name
WesternDiet   AGORA Western Diet   glc__D      D-Glucose
WesternDiet   AGORA Western Diet   ala__L      L-Alanine
WesternDiet   AGORA Western Diet   cpd00027    D-Glucose (SEED)
LB            Lysogeny broth       ala__L      L-Alanine
```

BiGG (`glc__D`) and ModelSEED (`cpd00027`) IDs are both accepted. With `--medium-name`, fast-mic translates BiGG → SEED via `--compounds-tsv` (ModelSEED `compounds.tsv`) so AGORA-style and gapseq-style models are both matched.

### CSV medium file / CSV 培养基文件

For gapseq models, use `--medium-file` with CSV format:

```
compounds,name,maxFlux
cpd00027,D-Glucose,10.0
cpd00035,L-Alanine,1.0
cpd00009,Phosphate,1000.0
```

### Tiered uptake limits / 分级摄取限制

Compound class is auto-detected from compound ID.

| Class | Default rate | Examples |
|---|---|---|
| Carbon sources | `--medium-uptake-limit` (10.0 mmol/gDW/h) | Sugars, organic acids |
| Amino acids | 1.0 mmol/gDW/h | Standard 20 + D-forms + ornithine |
| Nucleobases / nucleosides | 0.5 mmol/gDW/h | Adenine, uracil, adenosine |
| Cofactors / vitamins | 0.1 mmol/gDW/h | Folate, B12, riboflavin |
| Inorganic ions | Unlimited (−1000) | Na⁺, K⁺, Fe²⁺, PO₄³⁻ |

---

## Examples / 示例

### Pairwise + SCFA tracking / 两两配对 + 短链脂肪酸追踪

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

## Downstream analysis / 下游分析流水线

fast-mic produces TSV outputs that feed a set of post-processing scripts in `scripts/`. All scripts run from the repository root.

### 1. Prebiotic gradient / 益生元梯度分析

`media/` contains 10 gapseq medium CSV files (`gradient_L0_base_gapseq.csv` through `gradient_L9_mos_gapseq.csv`) spanning a cumulative prebiotic gradient. Each level adds one prebiotic's hydrolysis products on top of the mucin-containing base (design: Akkermansia viable at all levels).

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

**Run gradient** / 运行梯度分析:

```bash
bash scripts/run_gradient.sh \
     --group1 test/akk/akk_genomes_faa_gapseq_wdm_xml \
     --group2 test/UHGG/final_bacteria_gapseq_xml \
     --threads 12 --full-tsv
# Output: gradient_result/L{0-9}_*.tsv  +  L{0-9}_*.full.tsv
```

### 2. Paper figures & supplementary tables / 论文图表生成

**Generate paper figures** / 生成论文图（输出至 `results/figures_paper/`，TIFF/PDF/PNG）:

```bash
Rscript scripts/plot_fig1_benchmark.R                # Figure 1: fast-mic vs COBRApy benchmark (accuracy + thread scaling)
Rscript scripts/plot_fig2_phylo_trees.R              # Figure 2: phylogenetic trees (Akkermansia / Lactobacillus-group)
Rscript scripts/plot_fig3_gradient_overview.R        # Figure 3: interaction fractions across the L0–L9 gradient
Rscript scripts/plot_fig4_mechanisms.R               # Figure 4: competition intensity, cross-feeding, growth ratios
Rscript scripts/plot_fig5_crossfeed_sankey.R         # Figure 5: genome-encoded cross-feeding currencies (Sankey)
Rscript scripts/plot_fig6_strain_heatmap.R           # Figure 6: strain-resolved mutualism heatmaps
Rscript scripts/plot_figS1_interaction_composition.R # Fig S1: full six-category interaction composition
Rscript scripts/plot_figS2_benefit_landscape.R       # Fig S2: β_A vs β_B benefit landscape
Rscript scripts/plot_figS3_crossfeed_landscape.R     # Fig S3: top cross-fed metabolites per system
Rscript scripts/plot_figS4_threshold_sensitivity.R   # Fig S4: robustness to classification thresholds
Rscript scripts/plot_figS5_objective_sensitivity.R   # Fig S5: objective-function sensitivity (lex vs fixed-ratio)
```

**Generate supplementary tables** / 生成附表:

```bash
python3 scripts/make_supplementary_tables.py
# Output: results/figures_paper/Supplementary_Tables.xlsx
#   (genome panel & CheckM2 QC, prebiotic levels, gradient interaction summary,
#    strain-resolved mutualism, cross-fed metabolites, benchmark accuracy & scaling)
```

### 3. Archived analyses / 已归档分析

Earlier multi-site candidate analysis, EIR/EcoGS abundance-weighting, and HPLC/GC validation-hypothesis scripts have been moved to `scripts/_archive/` and are no longer part of the main pipeline. / 早期的多生境候选分析、EIR/EcoGS 丰度加权与 HPLC/GC 验证假说脚本已移至 `scripts/_archive/`，不再属于主流程。

---

## Benchmarking / 基准测试

A standalone binary `bench-single-fba` measures single-species FBA throughput across many models — under one medium or several at once — and supports thread-scaling benchmarks. / 独立二进制 `bench-single-fba` 在多个模型上测量单物种 FBA 通量，支持单个或多个培养基，以及线程扩展基准。

**Single medium / 单培养基** — a named medium from a TSV database, or one gapseq/SEED-format CSV:

```bash
# Named medium from a TSV database / TSV 库中的命名培养基
bench-single-fba media/media_db.tsv WesternDiet \
  --model-list models.txt --threads 0 > bench_results.tsv

# A single gapseq/SEED-format CSV medium / 单个 gapseq/SEED 格式 CSV 培养基
bench-single-fba --medium-file media/gradient_L0_base_gapseq.csv \
  --model-list models.txt --threads 0 > bench_results.tsv
```

**Several CSV media at once / 一次输入多个 CSV 培养基** — pass `--media-list FILE`, a plain-text file with one medium-CSV path per line (absolute or relative). Each model is loaded **once** and evaluated under **every** medium, yielding one row per (model, medium) pair. This is the workload behind the 10-level prebiotic gradient. / 通过 `--media-list FILE` 传入一个纯文本文件（每行一个培养基 CSV 路径，绝对或相对）。每个模型只加载一次，并在所有培养基下评估，输出每个（模型, 培养基）组合一行——这正是 10 级益生元梯度的工作负载。

```bash
# media/gradient_media_list.txt lists the 10 gradient CSVs (L0–L9), one path per line
bench-single-fba --media-list media/gradient_media_list.txt \
  --model-list models.txt --threads 0 > bench_gradient.tsv
```

Output columns: `model_id`, `n_metabolites`, `n_reactions`, `n_genes`, `biomass_rxn`, `growth_rate`, `load_time_s`, `fba_time_s`. With `--media-list`, `model_id` is suffixed with the medium label (the CSV file stem), e.g. `L_acidophilus_NCFM__gradient_L0_base_gapseq`, so each row stays uniquely keyed per (model, medium). / 输出列同上；使用 `--media-list` 时 `model_id` 会附加培养基标签（CSV 文件名主干），如 `L_acidophilus_NCFM__gradient_L0_base_gapseq`，保证每个（模型, 培养基）行唯一可索引。

### Reference validation against COBRApy / 与 COBRApy 对照验证

fast-mic ships an end-to-end correctness pipeline: `benchmark/run_thread_scaling.sh` runs both fast-mic and COBRApy (HiGHS backend) on the same model corpus and produces `benchmark_results/correctness/stats.tsv` with Pearson r, R², and MAE. A `cargo` test parses that file and asserts thresholds.

```bash
# 1. Regenerate the comparison TSVs (~hours on 1000 models; requires Python + cobrapy + highspy + osqp)
#    Single medium by default; set MEDIA_LIST to reproduce the 10-level gradient corpus above.
MEDIA_LIST=media/gradient_media_list.txt bash benchmark/run_thread_scaling.sh

# 2. Assert the contract: Pearson r ≥ 0.999, R² ≥ 0.998, MAE ≤ 1e-3
cargo test --test reference_validation -- --ignored
```

**Current status:** 9,950 genome–medium pairs (1,000 UHGG models across the L0–L9 gradient; 5 COBRApy timeouts excluded), Pearson r = 1.000, MAE = 3.12 × 10⁻⁷. Run this gate before tagging a release or merging changes that touch the LP path (`cobra.rs`, `medium.rs`, `sbml.rs`).

---

## Testing / 测试

```bash
# Default unit + integration tests (fast)
cargo test

# Reference validation (requires stats.tsv produced by benchmark above)
cargo test --test reference_validation -- --ignored

# Lint and format
cargo clippy -- -D warnings
cargo fmt -- --check
```

Test coverage:

| Module | Tests | Covers |
|---|---|---|
| `cobra` | 22 | Exchange detection, compartment inference, biomass finder, merged-model construction, interaction classification, NaN-safe sorts |
| `medium` | 4 | Cofactor pre-opened uptake preservation; non-cofactor closed; cofactor in medium uses tier bound; cofactor not pre-opened stays closed |
| Integration (`tests/`) | 1 (ignored) | COBRApy reference-agreement thresholds |

---

## Library API / 库 API

`fast-mic` is also usable as a Rust library.

```rust
use fast_mic::{cobra, medium, sbml};

let model = sbml::parse_sbml("model.xml")?;
let medium_set = medium::expand_medium_compounds(&base_compounds);
let params = cobra::AnalysisParams::default();
let result = cobra::run_fba(&model, &medium_set, &params)?;
println!("growth = {}", result.objective_value);
```

| Module | Purpose |
|---|---|
| `model` | Core types: `MetabolicModel`, `Reaction`, `Metabolite` (with optional `formula`), `PairwiseResult`, `InteractionType`. Also exports the unified `KNOWN_EXTERNAL_COMPARTMENTS`, `EXCHANGE_EXCLUDES`, `is_canonical_exchange`, `is_extracellular_compartment`. |
| `sbml` | SBML parsing (FBC v2). Two-pass, single read. Extracts `fbc:chemicalFormula`. |
| `medium` | Compound matching, exchange-reaction detection (COBRApy-style), tiered uptake bounds, cofactor pre-opened preservation, BiGG → SEED translation. |
| `cobra` | FBA, CycleFreeFlux + pFBA (`run_fba`, `run_fba_locked`), pairwise co-culture (lexicographic max-min, or fixed-ratio via `AnalysisParams::fixed_ratio`), cross-feeding analysis. `AnalysisParams::lock_tol` and `::fixed_ratio` configurable. |

---

## Citing / 引用

If you use fast-mic in published work, please cite:

- **Desouki *et al.* (2015)**, *CycleFreeFlux: efficient removal of thermodynamically infeasible loops from flux distributions*, BMC Bioinformatics **16**:283 — *cycle-removal algorithm.*
- **Bertsimas, Farias & Trichakis (2011)**, *The Price of Fairness*, Operations Research **59**(1):17-31 — *lexicographic max-min co-culture formulation.*

---

## License / 许可证

MIT or Apache-2.0, at your option.

## Issues & contributions / 问题反馈与贡献

Bug reports and PRs welcome at the project repository.
