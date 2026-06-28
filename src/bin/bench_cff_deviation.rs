//! Validation benchmark: post-CycleFreeFlux biomass deviation.
//!
//! For every (model, medium) the Stage-1 FBA optimum is compared with the
//! biomass returned after the CycleFreeFlux + pFBA (loop-removal) step. The
//! deviation `fba_optimal - post_cff_biomass` quantifies how far loop removal
//! moves the growth rate from the FBA optimum; it should not exceed the LP
//! solver tolerance. Reports the maximum deviation across all evaluations.
//!
//! Usage:
//!   bench-cff-deviation --media-list <file> --model-list <file> \
//!       [--threads N] [--medium-uptake-limit X]

use anyhow::Result;
use fast_mic::{cobra, medium, sbml};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;

struct Args {
    media_list: String,
    model_list: String,
    threads: usize,
    medium_uptake_limit: f64,
}

fn parse_args() -> Result<Args> {
    let raw: Vec<String> = std::env::args().collect();
    let (mut media_list, mut model_list) = (None, None);
    let (mut threads, mut lim) = (0usize, 10.0f64);
    let mut i = 1;
    while i < raw.len() {
        match raw[i].as_str() {
            "--media-list" => { i += 1; media_list = Some(raw[i].clone()); }
            "--model-list" => { i += 1; model_list = Some(raw[i].clone()); }
            "--threads" => { i += 1; threads = raw[i].parse()?; }
            "--medium-uptake-limit" => { i += 1; lim = raw[i].parse()?; }
            other => anyhow::bail!("unknown arg: {}", other),
        }
        i += 1;
    }
    Ok(Args {
        media_list: media_list.ok_or_else(|| anyhow::anyhow!("--media-list required"))?,
        model_list: model_list.ok_or_else(|| anyhow::anyhow!("--model-list required"))?,
        threads,
        medium_uptake_limit: lim,
    })
}

type Medium = (String, HashSet<String>, HashSet<String>, HashMap<String, f64>);

fn read_lines(path: &str) -> Result<Vec<String>> {
    Ok(std::fs::read_to_string(path)?
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .collect())
}

fn main() -> Result<()> {
    let args = parse_args()?;
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().ok();
    }

    let media: Vec<Medium> = read_lines(&args.media_list)?
        .iter()
        .map(|p| {
            let (compounds, per) = medium::parse_medium_csv(p)?;
            let label = std::path::Path::new(p)
                .file_stem().and_then(|s| s.to_str()).unwrap_or("medium").to_string();
            let expanded = medium::expand_medium_compounds(&compounds);
            Ok((label, compounds, expanded, per))
        })
        .collect::<Result<Vec<_>>>()?;

    let model_paths = read_lines(&args.model_list)?;
    eprintln!("Models: {} × {} media = {} evaluations",
        model_paths.len(), media.len(), model_paths.len() * media.len());

    let mub = medium::MediumBounds { carbon_source: args.medium_uptake_limit, ..Default::default() };
    let params = cobra::AnalysisParams::default();
    const MIN_VIABLE: f64 = 1e-4;

    // (model_id, medium_label, fba_optimal, post_cff_biomass, deviation, viable)
    let rows: Vec<(String, String, f64, f64, f64, bool)> = model_paths
        .par_iter()
        .flat_map(|path| {
            let mut out = Vec::new();
            if let Ok(mut model) = sbml::parse_sbml(path) {
                let id = model.id.clone();
                for (label, base, expanded, per) in &media {
                    medium::apply_medium(&mut model, base, &mub, Some(per));
                    if let Ok(r) = cobra::run_fba(&model, expanded, &params) {
                        let dev = r.fba_optimal - r.objective_value;
                        out.push((id.clone(), label.clone(), r.fba_optimal,
                                  r.objective_value, dev, r.objective_value >= MIN_VIABLE));
                    }
                }
            }
            out
        })
        .collect();

    // ── Per-row TSV (stdout) ──
    let stdout = std::io::stdout();
    let mut o = stdout.lock();
    writeln!(o, "model_id\tmedium\tfba_optimal\tpost_cff_biomass\tdeviation\tviable")?;
    for (m, med, fo, pc, d, v) in &rows {
        writeln!(o, "{}\t{}\t{:.10e}\t{:.10e}\t{:.10e}\t{}", m, med, fo, pc, d, v)?;
    }

    // ── Summary (stderr): max deviation over the GROWING (viable) evaluations ──
    let viable: Vec<f64> = rows.iter().filter(|r| r.5).map(|r| r.4.abs()).collect();
    let n_viable = viable.len();
    let max_dev = viable.iter().cloned().fold(0.0f64, f64::max);
    let mean_dev = if n_viable > 0 { viable.iter().sum::<f64>() / n_viable as f64 } else { 0.0 };
    let n_above_1em6 = viable.iter().filter(|&&d| d > 1e-6).count();
    // locate the worst case
    let worst = rows.iter().filter(|r| r.5)
        .max_by(|a, b| a.4.abs().partial_cmp(&b.4.abs()).unwrap());

    eprintln!();
    eprintln!("======================================================");
    eprintln!("  Post-CycleFreeFlux biomass deviation (validation)");
    eprintln!("======================================================");
    eprintln!("Total evaluations:        {}", rows.len());
    eprintln!("Viable (growing) evals:   {}", n_viable);
    eprintln!("HiGHS feasibility tol:    1e-7   |   lock_tol (ε): {:e}", params.lock_tol);
    eprintln!("MAX |deviation|:          {:.6e} h⁻¹", max_dev);
    eprintln!("Mean |deviation|:         {:.6e} h⁻¹", mean_dev);
    eprintln!("Evals with |dev| > 1e-6:  {}", n_above_1em6);
    if let Some(w) = worst {
        eprintln!("Worst case: {} @ {}  (FBA {:.8} → CFF {:.8}, Δ {:.3e})",
            w.0, w.1, w.2, w.3, w.4);
    }
    Ok(())
}
