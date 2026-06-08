use anyhow::Result;
use fast_mic::{cobra, medium, sbml};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

// ============================================================
// Benchmark result for one model
// ============================================================

struct ModelBenchResult {
    model_id: String,
    n_metabolites: usize,
    n_reactions: usize,
    n_genes: usize,
    biomass_rxn: String,
    growth_rate: f64,
    load_time_s: f64,
    fba_time_s: f64,
}

// ============================================================
// Media DB parsing
// ============================================================

fn parse_media_db(path: &str) -> Result<HashMap<String, (String, HashSet<String>)>> {
    let content = std::fs::read_to_string(path)
        .map_err(|e| anyhow::anyhow!("Cannot read media_db '{}': {}", path, e))?;
    let mut db: HashMap<String, (String, HashSet<String>)> = HashMap::new();
    for (i, line) in content.lines().enumerate() {
        if i == 0 {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let entry = db
            .entry(fields[0].trim().to_string())
            .or_insert_with(|| (fields[1].trim().to_string(), HashSet::new()));
        entry.1.insert(fields[2].trim().to_string());
    }
    Ok(db)
}

// ============================================================
// Peak memory (cross-platform) — returns KB internally
// ============================================================

#[cfg(target_os = "linux")]
fn get_peak_memory_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmPeak:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse().ok())
        })
        .unwrap_or(0)
}

#[cfg(target_os = "linux")]
fn get_rss_kb() -> u64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmRSS:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse().ok())
        })
        .unwrap_or(0)
}

#[cfg(target_os = "macos")]
fn get_peak_memory_kb() -> u64 {
    unsafe {
        let mut rusage: libc::rusage = std::mem::zeroed();
        if libc::getrusage(libc::RUSAGE_SELF, &mut rusage) == 0 {
            (rusage.ru_maxrss as u64) / 1024
        } else {
            0
        }
    }
}

#[cfg(target_os = "macos")]
fn get_rss_kb() -> u64 {
    get_peak_memory_kb()
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
fn get_peak_memory_kb() -> u64 {
    0
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
fn get_rss_kb() -> u64 {
    0
}

// ============================================================
// CLI argument parsing (manual — keeps binary self-contained)
// ============================================================

struct Args {
    media_db: String,
    medium_name: String,
    model_paths: Vec<String>,
    threads: usize,      // --threads N (0 = rayon default = all cores; 1 = serial)
    medium_file: Option<String>,   // --medium-file CSV (gapseq SEED-format)
    compounds_tsv: Option<String>, // --compounds-tsv (ModelSEED BiGG↔SEED aliases)
    medium_uptake_limit: f64,      // --medium-uptake-limit (carbon-source cap, mmol/gDW/h)
}

fn parse_args() -> Result<Args> {
    let raw: Vec<String> = std::env::args().collect();

    // Collect positional args (those not starting with --)
    let mut positional: Vec<String> = Vec::new();
    let mut model_list_path: Option<String> = None;
    let mut threads: usize = 0;
    let mut medium_file: Option<String> = None;
    let mut compounds_tsv: Option<String> = None;
    let mut medium_uptake_limit: f64 = 10.0; // matches fast-mic main binary default

    let mut i = 1usize;
    while i < raw.len() {
        match raw[i].as_str() {
            "--model-list" => {
                i += 1;
                model_list_path = Some(raw.get(i)
                    .ok_or_else(|| anyhow::anyhow!("--model-list needs a file path"))?
                    .clone());
            }
            "--threads" => {
                i += 1;
                threads = raw.get(i)
                    .ok_or_else(|| anyhow::anyhow!("--threads needs a value"))?
                    .parse()?;
            }
            "--medium-file" => {
                i += 1;
                medium_file = Some(raw.get(i)
                    .ok_or_else(|| anyhow::anyhow!("--medium-file needs a CSV path"))?
                    .clone());
            }
            "--compounds-tsv" => {
                i += 1;
                compounds_tsv = Some(raw.get(i)
                    .ok_or_else(|| anyhow::anyhow!("--compounds-tsv needs a path"))?
                    .clone());
            }
            "--medium-uptake-limit" => {
                i += 1;
                medium_uptake_limit = raw.get(i)
                    .ok_or_else(|| anyhow::anyhow!("--medium-uptake-limit needs a value"))?
                    .parse()?;
            }
            other if other.starts_with("--") => {
                anyhow::bail!("Unknown flag: {}", other);
            }
            _ => {
                positional.push(raw[i].clone());
            }
        }
        i += 1;
    }

    // When --medium-file is given, positional <media_db> <medium_name> are
    // optional (the CSV alone defines the medium). Otherwise they are required.
    let need_positional = if medium_file.is_some() { 0 } else { 2 };
    if positional.len() < need_positional {
        eprintln!("Usage: bench-single-fba <media_db.tsv> <medium_name> <model1.xml> [model2.xml ...]");
        eprintln!("       bench-single-fba <media_db.tsv> <medium_name> --model-list <file.txt>");
        eprintln!("       bench-single-fba --medium-file <medium.csv> --model-list <file.txt>");
        eprintln!("       [--threads N] [--compounds-tsv <compounds.tsv>]");
        eprintln!();
        eprintln!("  --model-list FILE     Read model paths from FILE (one per line), avoids ARG_MAX");
        eprintln!("  --medium-file CSV     Use a gapseq/SEED-format medium CSV (compounds,name,maxFlux)");
        eprintln!("                        Compatible with ModelSEED models (e.g. gapseq output).");
        eprintln!("                        When given, replaces <media_db> <medium_name> positionals.");
        eprintln!("  --compounds-tsv FILE  ModelSEED compounds.tsv — adds BiGG↔SEED alias translation");
        eprintln!("                        for media_db.tsv entries (helps BiGG media match SEED models)");
        eprintln!("  --threads N           Parallel worker threads (0=all cores, 1=serial; default 0)");
        eprintln!("  --medium-uptake-limit X   Carbon-source uptake cap (mmol/gDW/h, default 10.0)");
        eprintln!("                            Matches the fast-mic main binary's flag of the same name.");
        eprintln!();
        eprintln!("  Tiered uptake limits:");
        eprintln!("    carbon-source X (--medium-uptake-limit, default 10.0)");
        eprintln!("    amino-acid 1.0 / nucleobase 0.5 / cofactor 0.1 mmol/gDW/h");
        eprintln!("    (per-compound maxFlux from --medium-file CSV overrides these defaults)");
        std::process::exit(1);
    }

    // Build model path list: from --model-list file or from positional args
    let model_paths = if let Some(list_path) = model_list_path {
        let content = std::fs::read_to_string(&list_path)
            .map_err(|e| anyhow::anyhow!("Cannot read model list '{}': {}", list_path, e))?;
        let paths: Vec<String> = content.lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .collect();
        eprintln!("Read {} model paths from {}", paths.len(), list_path);
        paths
    } else {
        let model_start = need_positional;
        if positional.len() <= model_start {
            eprintln!("ERROR: Need at least one model path, or use --model-list <file>");
            std::process::exit(1);
        }
        positional[model_start..].to_vec()
    };

    // media_db / medium_name only required when --medium-file is NOT used
    let media_db = if need_positional >= 1 { positional[0].clone() } else { String::new() };
    let medium_name = if need_positional >= 2 { positional[1].clone() } else { String::new() };

    Ok(Args {
        media_db,
        medium_name,
        model_paths,
        threads,
        medium_file,
        compounds_tsv,
        medium_uptake_limit,
    })
}

// ============================================================
// Main
// ============================================================

fn main() -> Result<()> {
    let args = parse_args()?;

    // ── Resolve medium ──
    //
    // Two mutually-exclusive sources:
    //   (A) --medium-file <CSV>           — gapseq SEED-format CSV
    //                                        (compounds,name,maxFlux per row)
    //   (B) <media_db.tsv> <medium_name>  — legacy TSV database with BiGG IDs
    //
    // For (B), an optional --compounds-tsv can be supplied to translate BiGG
    // names to ModelSEED cpd IDs, so the same media_db.tsv works on both
    // BiGG (CarveMe/AGORA) and SEED (gapseq) models.
    let (medium_label, base_compounds, per_compound_bounds): (
        String,
        std::collections::HashSet<String>,
        Option<std::collections::HashMap<String, f64>>,
    ) = if let Some(ref csv_path) = args.medium_file {
        let (compounds, per_bounds) = medium::parse_medium_csv(csv_path)?;
        let label = std::path::Path::new(csv_path)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("medium-file")
            .to_string();
        (label, compounds, Some(per_bounds))
    } else {
        let media_db = parse_media_db(&args.media_db)?;
        let (desc, base_compounds) = media_db
            .get(&args.medium_name)
            .ok_or_else(|| {
                let available: Vec<_> = media_db.keys().collect();
                anyhow::anyhow!(
                    "Medium '{}' not found in '{}'. Available: {:?}",
                    args.medium_name,
                    args.media_db,
                    available
                )
            })?;

        // Optionally augment with BiGG→SEED aliases so SEED-format models
        // (e.g. gapseq) recognise BiGG-named medium entries.
        let mut compounds = base_compounds.clone();
        if let Some(ref tsv) = args.compounds_tsv {
            let map = medium::load_bigg_to_seed(tsv);
            let mut added = 0usize;
            for c in base_compounds.iter() {
                let stripped = c.strip_prefix("M_").unwrap_or(c);
                let root = stripped
                    .strip_suffix("_e0").or_else(|| stripped.strip_suffix("_e"))
                    .unwrap_or(stripped);
                if let Some(seeds) = map.get(root) {
                    for s in seeds {
                        if compounds.insert(format!("{}_e0", s)) { added += 1; }
                        compounds.insert(format!("{}_e", s));
                        compounds.insert(s.clone());
                    }
                }
            }
            eprintln!("  Compounds DB: {} → +{} SEED aliases from {}",
                map.len(), added, tsv);
        }
        (format!("{} ({})", args.medium_name, desc), compounds, None)
    };

    let medium_set = medium::expand_medium_compounds(&base_compounds);
    // Tiered uptake bounds (default: carbon 10 / aa 1 / nuc 0.5 / cofactor 0.1).
    // `--medium-uptake-limit` overrides the carbon-source cap, matching the
    // main fast-mic binary's flag of the same name.  Per-compound maxFlux
    // from a `--medium-file` CSV still overrides these defaults.
    let mub = medium::MediumBounds {
        carbon_source: args.medium_uptake_limit,
        ..Default::default()
    };

    eprintln!(
        "Medium: {} — {} base → {} expanded compounds",
        medium_label,
        base_compounds.len(),
        medium_set.len()
    );
    eprintln!(
        "Uptake limits: carbon={:.2}  aa={:.2}  nucleobase={:.2}  cofactor={:.2}{}",
        mub.carbon_source, mub.amino_acid, mub.nucleobase, mub.cofactor,
        if per_compound_bounds.is_some() {
            "  (per-compound overrides from CSV)"
        } else { "" }
    );
    eprintln!("Models: {}", args.model_paths.len());

    // ── Configure rayon thread pool ──
    // threads=0 means rayon default (num_cpus); 1 = serial; N = N workers.
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .map_err(|e| anyhow::anyhow!("Failed to init rayon pool: {}", e))?;
    }
    let effective_threads = rayon::current_num_threads();
    eprintln!("Threads: {} ({})",
        effective_threads,
        if args.threads == 0 { "default" } else { "user-specified" }
    );
    eprintln!();

    let params = cobra::AnalysisParams::default();

    let total_start = Instant::now();
    let rss_before = get_rss_kb();

    // ── Progress bar (parallel-safe via indicatif) ──
    let total = args.model_paths.len() as u64;
    let pb = ProgressBar::new(total);
    pb.set_style(
        ProgressStyle::with_template(
            "  [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({per_sec}, ETA {eta})"
        )
        .unwrap()
        .progress_chars("=>-"),
    );

    let counter = AtomicUsize::new(0);

    // ── Parallel per-model pFBA. par_iter preserves input order in collect(). ──
    let processed: Vec<(usize, Result<ModelBenchResult>)> = args.model_paths
        .par_iter()
        .enumerate()
        .map(|(idx, path)| {
            let res = (|| -> Result<ModelBenchResult> {
                let load_start = Instant::now();
                let mut model = sbml::parse_sbml(path)?;
                medium::apply_medium(&mut model, &base_compounds, &mub, per_compound_bounds.as_ref());
                let load_elapsed = load_start.elapsed();

                let fba_start = Instant::now();
                let fba_result = cobra::run_fba(&model, &medium_set, &params)?;
                let fba_elapsed = fba_start.elapsed();

                let biomass_id = model
                    .biomass_reaction
                    .as_deref()
                    .unwrap_or("NONE")
                    .to_string();

                Ok(ModelBenchResult {
                    model_id: model.id.clone(),
                    n_metabolites: model.metabolites.len(),
                    n_reactions: model.reactions.len(),
                    n_genes: model.genes.len(),
                    biomass_rxn: biomass_id,
                    growth_rate: fba_result.objective_value,
                    load_time_s: load_elapsed.as_secs_f64(),
                    fba_time_s: fba_elapsed.as_secs_f64(),
                })
            })();

            let n = counter.fetch_add(1, Ordering::Relaxed) + 1;
            pb.set_position(n as u64);
            if let Err(ref e) = res {
                pb.println(format!("  [{}/{}] FAILED {}: {}", n, total, path, e));
            }
            (idx, res)
        })
        .collect();

    pb.finish_and_clear();

    // ── Separate successes from failures, preserve original order ──
    let mut sorted = processed;
    sorted.sort_by_key(|(idx, _)| *idx);
    let mut results: Vec<ModelBenchResult> = Vec::with_capacity(sorted.len());
    let mut n_failed = 0usize;
    for (_, r) in sorted {
        match r {
            Ok(br) => results.push(br),
            Err(_) => n_failed += 1,
        }
    }
    if n_failed > 0 {
        eprintln!("  ⚠ {} models failed (see errors above)", n_failed);
    }

    let total_elapsed = total_start.elapsed();
    let peak_mem_kb = get_peak_memory_kb();
    let rss_after = get_rss_kb();

    // ── TSV output to stdout ──
    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    writeln!(
        out,
        "model_id\tn_metabolites\tn_reactions\tn_genes\tbiomass_rxn\tgrowth_rate\tload_time_s\tfba_time_s"
    )?;
    for r in &results {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}",
            r.model_id,
            r.n_metabolites,
            r.n_reactions,
            r.n_genes,
            r.biomass_rxn,
            r.growth_rate,
            r.load_time_s,
            r.fba_time_s,
        )?;
    }

    // ── Summary to stderr ──
    let n = results.len() as f64;
    let total_load_s: f64 = results.iter().map(|r| r.load_time_s).sum();
    let total_fba_s: f64 = results.iter().map(|r| r.fba_time_s).sum();
    let wall_s = total_elapsed.as_secs_f64();
    let peak_mem_mb = peak_mem_kb as f64 / 1024.0;
    let rss_delta_kb = rss_after.saturating_sub(rss_before);
    let rss_delta_mb = rss_delta_kb as f64 / 1024.0;

    eprintln!();
    eprintln!("========================================");
    eprintln!("  fast-mic Single-Species pFBA Benchmark");
    eprintln!("========================================");
    eprintln!("Models:              {}", results.len());
    eprintln!(
        "Total load time:     {:.4} s  ({:.4} s/model)",
        total_load_s,
        total_load_s / n
    );
    eprintln!(
        "Total pFBA time:     {:.4} s  ({:.4} s/model)",
        total_fba_s,
        total_fba_s / n
    );
    eprintln!(
        "Total wall time:     {:.4} s  ({:.4} s/model)",
        wall_s,
        wall_s / n
    );
    eprintln!("Peak memory:         {:.1} MB", peak_mem_mb);
    eprintln!("RSS delta:           {:.1} MB", rss_delta_mb);

    Ok(())
}
