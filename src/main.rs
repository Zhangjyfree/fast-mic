use anyhow::Result;
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;
use std::time::Instant;

mod model;
mod sbml;
mod medium;
mod cobra;

use model::*;

#[derive(Parser, Debug)]
#[command(name = "fast-mic", version, about, long_about = None)]
struct Cli {
    // ── Input ────────────────────────────────────────────────────────────
    #[arg(help_heading = "Input")]
    files: Vec<String>,

    #[arg(long, help_heading = "Input")]
    group1: Option<String>,

    #[arg(long, help_heading = "Input")]
    group2: Option<String>,

    /// Named medium from --media-db (e.g. "M9", "LB", "WesternDiet").
    /// Mutually exclusive with --medium-file.
    #[arg(long, help_heading = "Input")]
    medium_name: Option<String>,

    /// CSV medium file with SEED compound IDs and per-compound maxFlux
    /// (columns: compounds,name,maxFlux).  Mutually exclusive with --medium-name.
    #[arg(long, help_heading = "Input")]
    medium_file: Option<String>,

    #[arg(long, default_value = "media_db.tsv", help_heading = "Input")]
    media_db: String,

    /// ModelSEED compounds.tsv — used to translate BiGG compound IDs in
    /// --media-db to SEED IDs so gapseq models are matched correctly.
    /// Ignored when --medium-file is used (that file already contains SEED IDs).
    #[arg(long, default_value = "compounds.tsv", help_heading = "Input")]
    compounds_tsv: String,

    #[arg(long, help_heading = "Input")]
    pair_filter: Option<String>,

    // ── Output ───────────────────────────────────────────────────────────
    #[arg(short, long, default_value = "output.tsv", help_heading = "Output")]
    output: String,

    #[arg(long, help_heading = "Output")]
    full_tsv: Option<String>,

    #[arg(long, help_heading = "Output")]
    json: Option<String>,

    #[arg(short, long, help_heading = "Output")]
    verbose: bool,

    #[arg(long, help_heading = "Output")]
    summary: bool,

    // ── Medium uptake limits ─────────────────────────────────────────────
    #[arg(long, default_value_t = 10.0, help_heading = "Medium uptake limits")]
    medium_uptake_limit: f64,

    // ── Target-reaction tracking ─────────────────────────────────────────
    #[arg(long, value_delimiter = ',', help_heading = "Target-reaction tracking")]
    target_reaction: Vec<String>,

    // ── LP tolerances ────────────────────────────────────────────────────
    /// Tolerance for pinning reaction fluxes in CFF / pFBA lock constraints
    /// (`v ∈ [v* − tol, v* + tol]`).  Default 1e-5 (= 100× HiGHS feasibility
    /// tolerance); drop to 1e-7 only for exact single-species reproduction.
    #[arg(long, default_value_t = 1e-5, help_heading = "LP tolerances")]
    lock_tol: f64,


    // ── Performance ──────────────────────────────────────────────────────
    #[arg(long, default_value_t = 0, help_heading = "Performance")]
    threads: usize,

    #[arg(long, default_value_t = true, help_heading = "Performance")]
    cache_monoculture: bool,

    /// Use the fixed-ratio co-culture objective (μ_A/μ_B pinned to the
    /// monoculture ratio) instead of the default lexicographic max-min
    /// allocation.  Intended for objective-sensitivity analysis.
    #[arg(long, help_heading = "Co-culture objective")]
    fixed_ratio: bool,

    // ── DEPRECATED (kept hidden for backwards CLI compatibility) ────────
    #[arg(long, hide = true)]
    max_synergy: Option<f64>,
    #[arg(long, hide = true)]
    max_exchange_ratio: Option<f64>,
    #[arg(long, hide = true)]
    geneless_cap: Option<f64>,
}

// ============================================================
// Media DB parsing
// ============================================================

fn parse_media_db(path: &str) -> Result<HashMap<String, (String, HashSet<String>)>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("Cannot read media_db '{}': {}", path, e))?;
    let mut db: HashMap<String, (String, HashSet<String>)> = HashMap::new();
    for (i, line) in content.lines().enumerate() {
        if i == 0 { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 { continue; }
        let entry = db.entry(fields[0].trim().to_string())
            .or_insert_with(|| (fields[1].trim().to_string(), HashSet::new()));
        entry.1.insert(fields[2].trim().to_string());
    }
    Ok(db)
}

/// Resolved medium: compound set, expanded set, and optional per-compound bounds.
struct MediumSpec {
    base_compounds:     HashSet<String>,
    expanded_compounds: HashSet<String>,
    /// Per-compound max uptake (mmol/gDW/h).  Present when loaded from a CSV
    /// file.  `None` when loaded from the TSV media database (tiered bounds
    /// are used instead).
    per_compound_bounds: Option<HashMap<String, f64>>,
}

fn build_medium(cli: &Cli) -> Result<MediumSpec> {
    match (&cli.medium_file, &cli.medium_name) {
        (Some(_), Some(_)) => anyhow::bail!(
            "--medium-file and --medium-name are mutually exclusive; provide only one."),
        (None, None) => anyhow::bail!(
            "Provide either --medium-file <path.csv> or --medium-name <name>."),

        // ── CSV medium file (SEED IDs + explicit maxFlux per compound) ─────
        (Some(path), None) => {
            let (compounds, per_bounds) = medium::parse_medium_csv(path)?;
            let expanded = medium::expand_medium_compounds(&compounds);
            let label = std::path::Path::new(path)
                .file_stem().map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| path.clone());
            println!("📋 Medium: {} (CSV) — {} compounds → {} expanded (per-compound maxFlux)",
                label, compounds.len(), expanded.len());
            Ok(MediumSpec {
                base_compounds: compounds,
                expanded_compounds: expanded,
                per_compound_bounds: Some(per_bounds),
            })
        }

        // ── TSV media database (BiGG IDs, tiered uptake bounds) ───────────
        (None, Some(name)) => {
            let media_db = parse_media_db(&cli.media_db)?;
            if let Some((desc, base)) = media_db.get(name) {
                let bigg_to_seed = medium::load_bigg_to_seed(&cli.compounds_tsv);
                let mut augmented = base.clone();
                let mut n_seed_added = 0usize;
                for bigg_id in base.iter() {
                    if let Some(seed_ids) = bigg_to_seed.get(bigg_id.as_str()) {
                        for seed_id in seed_ids {
                            if augmented.insert(seed_id.clone()) {
                                n_seed_added += 1;
                            }
                        }
                    }
                }
                let expanded = medium::expand_medium_compounds(&augmented);
                println!("📋 Medium: {} ({}) — {} base → {} expanded{}",
                    name, desc, base.len(), expanded.len(),
                    if n_seed_added > 0 { format!(" (+{} SEED IDs)", n_seed_added) } else { String::new() });
                Ok(MediumSpec {
                    base_compounds: augmented,
                    expanded_compounds: expanded,
                    per_compound_bounds: None,
                })
            } else {
                let available: Vec<_> = media_db.keys().collect();
                anyhow::bail!("Medium '{}' not found in '{}'. Available: {:?}",
                    name, cli.media_db, available);
            }
        }
    }
}

fn load_pair_filter(path: &str) -> Result<HashSet<(String, String)>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("Cannot read pair filter '{}': {}", path, e))?;
    let mut pairs = HashSet::new();
    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 2 { continue; }
        if fields[0].eq_ignore_ascii_case("species_a") { continue; }
        pairs.insert((fields[0].trim().to_string(), fields[1].trim().to_string()));
    }
    Ok(pairs)
}

fn load_models_from_files(
    paths: &[String], base_compounds: &HashSet<String>,
    medium_bounds: &medium::MediumBounds,
    per_compound_bounds: Option<&HashMap<String, f64>>,
    verbose: bool,
) -> Vec<MetabolicModel> {
    let pb = ProgressBar::new(paths.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{bar:40.cyan/blue}] {pos}/{len} {msg}").unwrap()
        .progress_chars("#>-"));

    let mut indexed: Vec<(usize, Option<MetabolicModel>)> = paths.par_iter().enumerate()
        .map(|(idx, path)| {
            let res = match sbml::parse_sbml(path) {
                Ok(mut model) => {
                    medium::apply_medium(&mut model, base_compounds, medium_bounds, per_compound_bounds);
                    if verbose {
                        pb.println(format!(
                            "  ✓ {}: {} metabolites, {} reactions, {} genes",
                            model.id, model.metabolites.len(),
                            model.reactions.len(), model.genes.len()));
                    }
                    Some(model)
                }
                Err(e) => {
                    pb.println(format!("  ✗ Failed to parse {}: {}", path, e));
                    None
                }
            };
            pb.inc(1);
            (idx, res)
        }).collect();

    indexed.sort_by_key(|(idx, _)| *idx);
    let models: Vec<MetabolicModel> = indexed.into_iter().filter_map(|(_, m)| m).collect();
    pb.finish_with_message("Done");
    models
}

fn load_models_from_dir(
    dir: &str, base_compounds: &HashSet<String>,
    medium_bounds: &medium::MediumBounds,
    per_compound_bounds: Option<&HashMap<String, f64>>,
    verbose: bool,
) -> Result<Vec<MetabolicModel>> {
    let mut paths: Vec<String> = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map(|e| e == "xml" || e == "sbml").unwrap_or(false) {
            paths.push(path.to_string_lossy().to_string());
        }
    }
    paths.sort();
    if paths.is_empty() { anyhow::bail!("No .xml/.sbml files found in {}", dir); }
    println!("  📁 {} → {} files", dir, paths.len());
    Ok(load_models_from_files(&paths, base_compounds, medium_bounds, per_compound_bounds, verbose))
}

fn generate_pairs(n_group1: usize, n_group2: Option<usize>) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    match n_group2 {
        Some(n2) => for i in 0..n_group1 {
            for j in n_group1..(n_group1 + n2) { pairs.push((i, j)); }
        },
        None => for i in 0..n_group1 {
            for j in (i + 1)..n_group1 { pairs.push((i, j)); }
        },
    }
    pairs
}

fn precompute_monocultures(
    models: &[MetabolicModel], medium: &HashSet<String>, params: &cobra::AnalysisParams,
) -> HashMap<String, cobra::FBAResult> {
    let pb = ProgressBar::new(models.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} mono-pFBA")
        .unwrap().progress_chars("#>-"));

    let cache: Vec<(String, cobra::FBAResult)> = models.par_iter()
        .map(|m| {
            let res = match cobra::run_fba(m, medium, params) {
                Ok(r) => r,
                Err(e) => {
                    pb.println(format!("  ⚠ mono-pFBA failed for {}: {}", m.id, e));
                    cobra::FBAResult { objective_value: 0.0, flux_distribution: HashMap::new() }
                }
            };
            pb.inc(1);
            (m.id.clone(), res)
        }).collect();
    pb.finish_with_message("Done");
    cache.into_iter().collect()
}

fn fmt_frac(v: f64) -> String {
    if v.is_nan() { "NaN".to_string() } else { format!("{:.4}", v) }
}

fn write_tsv(path: &str, results: &[PairwiseResult], target_reactions: &[String]) -> Result<()> {
    use std::io::Write;
    let mut file = std::fs::File::create(path)?;

    let mut header: Vec<String> = [
        "species_a","species_b",
        "growth_a_alone","growth_b_alone",
        "growth_a_co","growth_b_co",
        "benefit_a","benefit_b",
        "interaction_type",
        "gene_supported_fraction",
        "n_exchanged_metabolites",
        "competition_intensity",
    ].iter().map(|s| s.to_string()).collect();

    for tr in target_reactions {
        header.push(format!("{}__alone_a", tr));
        header.push(format!("{}__alone_b", tr));
        header.push(format!("{}__co_total", tr));
        header.push(format!("{}__co_a", tr));
        header.push(format!("{}__co_b", tr));
    }
    writeln!(file, "{}", header.join("\t"))?;

    for r in results {
        let n_ex = r.a_to_b_exchanges.len() + r.b_to_a_exchanges.len();
        let mut row = vec![
            r.species_a.clone(), r.species_b.clone(),
            format!("{:.6}", r.growth_a_alone), format!("{:.6}", r.growth_b_alone),
            format!("{:.6}", r.growth_a_co), format!("{:.6}", r.growth_b_co),
            format!("{:.6}", r.benefit_a), format!("{:.6}", r.benefit_b),
            r.interaction_type.to_string(),
            fmt_frac(r.gene_supported_fraction),
            n_ex.to_string(),
            format!("{:.4}", r.competition_intensity),
        ];
        for tr in target_reactions {
            match r.target_fluxes.iter().find(|t| &t.reaction_id == tr) {
                Some(rec) => {
                    row.push(format!("{:.6}", rec.alone_a));
                    row.push(format!("{:.6}", rec.alone_b));
                    row.push(format!("{:.6}", rec.co_total));
                    row.push(format!("{:.6}", rec.co_a_contribution));
                    row.push(format!("{:.6}", rec.co_b_contribution));
                }
                None => for _ in 0..5 { row.push("NA".to_string()); }
            }
        }
        writeln!(file, "{}", row.join("\t"))?;
    }
    Ok(())
}

fn write_full_tsv(path: &str, results: &[PairwiseResult]) -> Result<()> {
    use std::io::Write;
    let mut file = std::fs::File::create(path)?;
    // ── Header ──
    // a_to_b_inferred / b_to_a_inferred are `;`-separated booleans
    // aligned 1:1 with a_to_b_metabolites / b_to_a_metabolites:
    //   "false" → direct cross-feeding (mass-balance evidence in the
    //              co-culture LP solution).
    //   "true"  → inferred / facilitation entry, attributed via the
    //              receiver's Δμ > 0 fallback (hypothesis-grade).
    let header = [
        "species_a","species_b",
        "growth_a_alone","growth_b_alone",
        "growth_a_co","growth_b_co",
        "mu","reciprocity_index",
        "benefit_a","benefit_b","net_interaction",
        "interaction_type","gene_supported_fraction",
        "a_to_b_metabolites","a_to_b_fluxes",
        "a_to_b_donor_genes","a_to_b_receiver_genes",
        "a_to_b_inferred","a_to_b_low_confidence",
        "b_to_a_metabolites","b_to_a_fluxes",
        "b_to_a_donor_genes","b_to_a_receiver_genes",
        "b_to_a_inferred","b_to_a_low_confidence",
        "competed_metabolites","competed_uptake_a","competed_uptake_b",
        "competed_genes_a","competed_genes_b",
        "competition_intensity",
    ];
    writeln!(file, "{}", header.join("\t"))?;
    for r in results {
        let a2b_mets = fmt_list(r.a_to_b_exchanges.iter().map(|e| &e.metabolite_id));
        let a2b_fluxes = r.a_to_b_exchanges.iter().map(|e| format!("{:.4}", e.flux))
            .collect::<Vec<_>>().join(";");
        let a2b_dgenes = r.a_to_b_exchanges.iter().map(|e| e.donor_genes.join(","))
            .collect::<Vec<_>>().join(";");
        let a2b_rgenes = r.a_to_b_exchanges.iter().map(|e| e.receiver_genes.join(","))
            .collect::<Vec<_>>().join(";");
        // Empty string (not "none") for empty lists -- see fmt_list().
        let a2b_infer = r.a_to_b_exchanges.iter().map(|e| e.inferred.to_string())
            .collect::<Vec<_>>().join(";");
        let a2b_lc = r.a_to_b_exchanges.iter().map(|e| e.low_confidence.to_string())
            .collect::<Vec<_>>().join(";");
        let b2a_mets = fmt_list(r.b_to_a_exchanges.iter().map(|e| &e.metabolite_id));
        let b2a_fluxes = r.b_to_a_exchanges.iter().map(|e| format!("{:.4}", e.flux))
            .collect::<Vec<_>>().join(";");
        let b2a_dgenes = r.b_to_a_exchanges.iter().map(|e| e.donor_genes.join(","))
            .collect::<Vec<_>>().join(";");
        let b2a_rgenes = r.b_to_a_exchanges.iter().map(|e| e.receiver_genes.join(","))
            .collect::<Vec<_>>().join(";");
        let b2a_infer = r.b_to_a_exchanges.iter().map(|e| e.inferred.to_string())
            .collect::<Vec<_>>().join(";");
        let b2a_lc = r.b_to_a_exchanges.iter().map(|e| e.low_confidence.to_string())
            .collect::<Vec<_>>().join(";");
        let comp_mets = fmt_list(r.shared_uptakes.iter().map(|s| &s.metabolite_id));
        let comp_uptake_a = r.shared_uptakes.iter().map(|s| format!("{:.4}", s.uptake_a))
            .collect::<Vec<_>>().join(";");
        let comp_uptake_b = r.shared_uptakes.iter().map(|s| format!("{:.4}", s.uptake_b))
            .collect::<Vec<_>>().join(";");
        let comp_ga = r.shared_uptakes.iter().map(|s| s.genes_a.join(","))
            .collect::<Vec<_>>().join(";");
        let comp_gb = r.shared_uptakes.iter().map(|s| s.genes_b.join(","))
            .collect::<Vec<_>>().join(";");
        let row = vec![
            r.species_a.clone(), r.species_b.clone(),
            format!("{:.6}", r.growth_a_alone), format!("{:.6}", r.growth_b_alone),
            format!("{:.6}", r.growth_a_co), format!("{:.6}", r.growth_b_co),
            format!("{:.4}", r.mu), format!("{:.4}", r.reciprocity_index),
            format!("{:.6}", r.benefit_a), format!("{:.6}", r.benefit_b),
            format!("{:.6}", r.net_interaction),
            r.interaction_type.to_string(),
            fmt_frac(r.gene_supported_fraction),
            a2b_mets, a2b_fluxes, a2b_dgenes, a2b_rgenes, a2b_infer, a2b_lc,
            b2a_mets, b2a_fluxes, b2a_dgenes, b2a_rgenes, b2a_infer, b2a_lc,
            comp_mets, comp_uptake_a, comp_uptake_b, comp_ga, comp_gb,
            format!("{:.4}", r.competition_intensity),
        ];
        writeln!(file, "{}", row.join("\t"))?;
    }
    Ok(())
}

/// Join a list of strings with `;` for TSV output.
///
/// Returns an empty string (NOT the literal "none") when the list is empty
/// so that downstream `split(";")` + `explode` pipelines don't manufacture
/// a fake `"none"` metabolite node. Empty fields are the canonical TSV
/// representation of "this pair had no exchanges in this direction".
fn fmt_list<'a, I, S>(iter: I) -> String
where I: Iterator<Item = &'a S>, S: AsRef<str> + 'a,
{
    let v: Vec<&str> = iter.map(|s| s.as_ref()).collect();
    if v.is_empty() { String::new() } else { v.join(";") }
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    println!("🦠 Fast Metabolic Interaction Calculator v{}", env!("CARGO_PKG_VERSION"));
    println!("================================================");
    println!("Mode: pairwise (FBA + CycleFreeFlux + lexicographic max-min)");

    if cli.max_synergy.is_some() || cli.max_exchange_ratio.is_some() || cli.geneless_cap.is_some() {
        eprintln!("  ⚠ --max-synergy / --max-exchange-ratio / --geneless-cap are deprecated and ignored.");
        eprintln!("    Empirical synergy/ratio caps were replaced by CycleFreeFlux + lex max-min.");
    }

    let start = Instant::now();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .map_err(|e| anyhow::anyhow!("Failed to init rayon pool: {}", e))?;
    }
    let effective_threads = rayon::current_num_threads();
    println!("⚙️  Threads: {} ({})",
        effective_threads,
        if cli.threads == 0 { "default" } else { "user-specified" });

    let params = cobra::AnalysisParams {
        lock_tol: cli.lock_tol,
        fixed_ratio: cli.fixed_ratio,
    };

    if !cli.summary {
        println!("⚙️  Medium uptake limits: carbon={} aa=1.0 mmol/gDW/h",
            cli.medium_uptake_limit);
        println!("⚙️  Monoculture pFBA cache: {}",
            if cli.cache_monoculture { "ENABLED" } else { "disabled" });
    }

    let medium_spec = build_medium(&cli)?;
    let base_compounds  = &medium_spec.base_compounds;
    let medium_compounds = &medium_spec.expanded_compounds;
    let per_cpd = medium_spec.per_compound_bounds.as_ref();
    let mub = medium::MediumBounds {
        carbon_source: cli.medium_uptake_limit,
        ..Default::default()
    };

    println!("\n📥 Loading models...");
    let (all_models, n_group1, n_group2) = if cli.group1.is_some() {
        let g1 = load_models_from_dir(
            cli.group1.as_ref().unwrap(), base_compounds, &mub, per_cpd, cli.verbose)?;
        let n1 = g1.len();
        if let Some(ref g2_dir) = cli.group2 {
            let g2 = load_models_from_dir(g2_dir, base_compounds, &mub, per_cpd, cli.verbose)?;
            let n2 = g2.len();
            let mut all = g1; all.extend(g2);
            (all, n1, Some(n2))
        } else { (g1, n1, None) }
    } else if !cli.files.is_empty() {
        let models = load_models_from_files(&cli.files, base_compounds, &mub, per_cpd, cli.verbose);
        let n = models.len();
        (models, n, None)
    } else {
        anyhow::bail!("Provide model files or use --group1 [--group2]");
    };

    if all_models.len() < 2 { anyhow::bail!("Need at least 2 models for pairwise analysis"); }

    if !cli.summary {
        println!("\n📊 Model Summary:");
        for m in &all_models {
            let bio = m.biomass_reaction.as_deref().unwrap_or("NONE");
            println!("  {}: {} metabolites, {} reactions, {} genes, biomass: {}",
                m.id, m.metabolites.len(), m.reactions.len(), m.genes.len(), bio);
        }
    }

    let mono_cache: HashMap<String, cobra::FBAResult> = if cli.cache_monoculture {
        println!("\n🧪 Pre-computing monoculture pFBA for all {} models...", all_models.len());
        let mono_start = Instant::now();
        let cache = precompute_monocultures(&all_models, &medium_compounds, &params);
        println!("  ✓ Done in {:.2}s ({} viable, {} non-viable)",
            mono_start.elapsed().as_secs_f64(),
            cache.values().filter(|r| r.objective_value > 1e-4).count(),
            cache.values().filter(|r| r.objective_value <= 1e-4).count());
        cache
    } else { HashMap::new() };

    let mut pairs = generate_pairs(n_group1, n_group2);
    let before_filter = pairs.len();

    if let Some(ref filter_path) = cli.pair_filter {
        let wanted = load_pair_filter(filter_path)?;
        let model_ids: Vec<&str> = all_models.iter().map(|m| m.id.as_str()).collect();
        pairs.retain(|&(i, j)| {
            let a = model_ids[i].to_string();
            let b = model_ids[j].to_string();
            wanted.contains(&(a.clone(), b.clone())) || wanted.contains(&(b, a))
        });
        println!("🔍 Pair filter '{}': {} → {} pairs", filter_path, before_filter, pairs.len());
        if pairs.is_empty() {
            anyhow::bail!("Pair filter matched 0 pairs — check model IDs in '{}'", filter_path);
        }
    }

    println!("\n🔬 Analyzing {} pairwise interactions...", pairs.len());

    let pb = ProgressBar::new(pairs.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} pairs")
        .unwrap().progress_chars("#>-"));
    let print_lock = Mutex::new(());

    let mut indexed: Vec<(usize, Result<PairwiseResult>)> = pairs.par_iter().enumerate()
        .map(|(idx, &(i, j))| {
            let id_a = all_models[i].id.as_str();
            let id_b = all_models[j].id.as_str();
            let cached_a = mono_cache.get(id_a);
            let cached_b = mono_cache.get(id_b);

            let res = cobra::calculate_pairwise_interaction(
                &all_models[i], &all_models[j],
                cached_a, cached_b,
                &medium_compounds, &params, &cli.target_reaction);

            if !cli.summary {
                if let Ok(ref result) = res {
                    let n_ex = result.a_to_b_exchanges.len() + result.b_to_a_exchanges.len();
                    let n_inferred = result.a_to_b_exchanges.iter().filter(|e| e.inferred).count()
                        + result.b_to_a_exchanges.iter().filter(|e| e.inferred).count();
                    let _g = print_lock.lock().unwrap();
                    let gsf_str = if result.gene_supported_fraction.is_nan() {
                        "NA".to_string()
                    } else {
                        format!("{:.2}", result.gene_supported_fraction)
                    };
                    pb.println(format!(
                        "  {} ↔ {} → {} (benefit: A={:+.1}%, B={:+.1}%) [{} ex ({} inferred), cf={:.2}, gsf={}]",
                        result.species_a, result.species_b,
                        result.interaction_type,
                        result.benefit_a * 100.0, result.benefit_b * 100.0,
                        n_ex, n_inferred, result.cross_feeding_score, gsf_str));
                    if cli.verbose {
                        pb.println(format!("    Single:  A={:.4}, B={:.4}",
                            result.growth_a_alone, result.growth_b_alone));
                        pb.println(format!("    Co-cul:  A={:.4}, B={:.4}",
                            result.growth_a_co, result.growth_b_co));
                        if !result.a_to_b_exchanges.is_empty() {
                            pb.println(format!("    {}→{}:", result.species_a, result.species_b));
                            for e in &result.a_to_b_exchanges {
                                let tag = if e.inferred { " [inferred]" } else { "" };
                                pb.println(format!("      {} (flux={:.3}){} [{}] → [{}]",
                                    e.metabolite_id, e.flux, tag,
                                    e.donor_genes.join(", "), e.receiver_genes.join(", ")));
                            }
                        }
                        if !result.b_to_a_exchanges.is_empty() {
                            pb.println(format!("    {}→{}:", result.species_b, result.species_a));
                            for e in &result.b_to_a_exchanges {
                                let tag = if e.inferred { " [inferred]" } else { "" };
                                pb.println(format!("      {} (flux={:.3}){} [{}] → [{}]",
                                    e.metabolite_id, e.flux, tag,
                                    e.donor_genes.join(", "), e.receiver_genes.join(", ")));
                            }
                        }
                        if !result.shared_uptakes.is_empty() {
                            pb.println("    Competed:".to_string());
                            for s in &result.shared_uptakes {
                                pb.println(format!("      {} (A:{:.2}/B:{:.2})",
                                    s.metabolite_id, s.uptake_a, s.uptake_b));
                            }
                        }
                    }
                }
            }
            pb.inc(1);
            (idx, res)
        }).collect();

    pb.finish_with_message("Done");
    indexed.sort_by_key(|(idx, _)| *idx);
    let mut results: Vec<PairwiseResult> = Vec::with_capacity(indexed.len());
    let mut n_failed = 0usize;
    for (_, r) in indexed {
        match r {
            Ok(pr) => results.push(pr),
            Err(e) => { eprintln!("  ⚠ Pair failed: {}", e); n_failed += 1; }
        }
    }
    if n_failed > 0 {
        eprintln!("  ⚠ {} pairs failed; {} succeeded.", n_failed, results.len());
    }

    write_tsv(&cli.output, &results, &cli.target_reaction)?;
    println!("\n💾 Results saved to {}", cli.output);
    if let Some(ref full_path) = cli.full_tsv {
        write_full_tsv(full_path, &results)?;
        println!("💾 Full TSV saved to {}", full_path);
    }
    if let Some(ref json_path) = cli.json {
        let json = serde_json::to_string_pretty(&results)?;
        std::fs::write(json_path, json)?;
        println!("💾 JSON saved to {}", json_path);
    }

    println!("\n📊 Summary:");
    let mut type_counts: HashMap<String, usize> = HashMap::new();
    for r in &results { *type_counts.entry(r.interaction_type.to_string()).or_insert(0) += 1; }
    let mut type_vec: Vec<_> = type_counts.iter().collect();
    type_vec.sort_by_key(|(_, c)| std::cmp::Reverse(**c));
    for (t, c) in &type_vec { println!("  {}: {} pairs", t, c); }

    if !results.is_empty() {
        let avg_a: f64 = results.iter().map(|r| r.benefit_a).sum::<f64>() / results.len() as f64;
        let avg_b: f64 = results.iter().map(|r| r.benefit_b).sum::<f64>() / results.len() as f64;

        // ── Direct vs. inferred cross-feeding tally ──
        let mut n_direct = 0usize;
        let mut n_inferred = 0usize;
        let mut flux_direct = 0.0f64;
        let mut flux_inferred = 0.0f64;
        for r in &results {
            for e in r.a_to_b_exchanges.iter().chain(r.b_to_a_exchanges.iter()) {
                if e.inferred { n_inferred += 1; flux_inferred += e.flux; }
                else { n_direct += 1; flux_direct += e.flux; }
            }
        }
        let n_total_ex = n_direct + n_inferred;
        let total_flux = flux_direct + flux_inferred;

        let valid_gsf: Vec<f64> = results.iter()
            .map(|r| r.gene_supported_fraction)
            .filter(|v| v.is_finite())
            .collect();
        let n_undefined = results.len() - valid_gsf.len();

        println!("  avg benefit_a: {:+.02}%, avg benefit_b: {:+.2}%", avg_a*100.0, avg_b*100.0);
        if n_total_ex > 0 {
            let pct_n = 100.0 * n_inferred as f64 / n_total_ex as f64;
            let pct_f = if total_flux > 0.0 { 100.0 * flux_inferred / total_flux } else { 0.0 };
            println!(
                "  cross-feeding entries: {} total ({} direct, {} inferred = {:.1}% of entries, {:.1}% of flux)",
                n_total_ex, n_direct, n_inferred, pct_n, pct_f,
            );
        } else {
            println!("  cross-feeding entries: 0 (no detected exchanges)");
        }
        if valid_gsf.is_empty() {
            println!("  gene_supported_fraction: NA (no pair had cross-feeding flux; n_undefined={})", n_undefined);
        } else {
            let mean_gsf: f64 = valid_gsf.iter().sum::<f64>() / valid_gsf.len() as f64;
            let n_full   = valid_gsf.iter().filter(|&&v| v >= 0.99).count();
            let n_high   = valid_gsf.iter().filter(|&&v| v >= 0.9 && v < 0.99).count();
            let n_mid    = valid_gsf.iter().filter(|&&v| v >= 0.1 && v < 0.9).count();
            let n_none   = valid_gsf.iter().filter(|&&v| v < 0.1).count();
            println!(
                "  gene_supported_fraction: mean={:.3} (n={}), \
                 ≥0.99: {}, [0.9,0.99): {}, [0.1,0.9): {}, <0.1: {}, NaN: {}",
                mean_gsf, valid_gsf.len(), n_full, n_high, n_mid, n_none, n_undefined,
            );
        }
    }

    // Report how many candidate exchange metabolites were filtered as
    // artifacts (cofactors, reactive intermediates, electron carriers, acyl-CoAs)
    // vs kept as real cross-feeding -- useful for methods section and review.
    cobra::print_filter_summary();

    let duration = start.elapsed();
    println!("\n⏱️  Total time: {:.2} s", duration.as_secs_f64());
    println!("✅ Done!");
    Ok(())
}
