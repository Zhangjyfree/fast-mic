// src/sbml.rs — IO-once two-pass version (recommended drop-in)
use crate::model::*;
use anyhow::{Context, Result};
use quick_xml::events::Event;
use quick_xml::Reader;
use std::collections::HashMap;
use std::path::Path;

struct GeneAssocFrame { operator: String, children: Vec<String> }

#[derive(Debug, Clone, Default)]
struct FbcObjective { reaction_coefficients: Vec<(String, f64)> }

pub fn parse_sbml<P: AsRef<Path>>(path: P) -> Result<MetabolicModel> {
    // ── Read file ONCE into memory; run two in-memory parses on it. ──
    let content = std::fs::read_to_string(&path)
        .with_context(|| format!("Failed to read {}", path.as_ref().display()))?;

    let mut model = MetabolicModel {
        id: path.as_ref().file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string(),
        name: None,
        metabolites: HashMap::new(),
        reactions: HashMap::new(),
        biomass_reaction: None,
        genes: HashMap::new(),
    };

    // ── Pass 1: parameters, gene products, FBC objectives ──
    let mut parameters: HashMap<String, f64> = HashMap::new();
    let mut active_objective_id: Option<String> = None;
    let mut objectives: HashMap<String, FbcObjective> = HashMap::new();
    let mut current_objective_id: Option<String> = None;

    {
        let mut xr = Reader::from_str(&content);
        xr.trim_text(true);
        let mut buf = Vec::new();
        loop {
            match xr.read_event_into(&mut buf) {
                Ok(Event::Empty(e)) | Ok(Event::Start(e)) => {
                    let tag_ref = e.name().as_ref().to_vec();
                    let tag_local = local_name(&tag_ref);
                    match tag_local {
                        b"parameter" => {
                            let mut id = String::new(); let mut value = 0.0;
                            for attr in e.attributes().flatten() {
                                match attr.key.as_ref() {
                                    b"id" => id = String::from_utf8_lossy(&attr.value).into_owned(),
                                    b"value" => {
                                        if let Ok(v) = std::str::from_utf8(&attr.value).unwrap_or("").parse::<f64>() {
                                            value = v;
                                        }
                                    }
                                    _ => {}
                                }
                            }
                            if !id.is_empty() { parameters.insert(id, value); }
                        }
                        b"geneProduct" => {
                            let (mut gid, mut gname, mut glabel) = (String::new(), String::new(), String::new());
                            for attr in e.attributes().flatten() {
                                let key = attr.key.as_ref();
                                let local = local_name(key);
                                let val = String::from_utf8_lossy(&attr.value).into_owned();
                                match local { b"id" => gid = val, b"name" => gname = val, b"label" => glabel = val, _ => {} }
                            }
                            if !gid.is_empty() {
                                model.genes.insert(gid.clone(), Gene { id: gid, name: gname, label: glabel });
                            }
                        }
                        b"listOfObjectives" => {
                            for attr in e.attributes().flatten() {
                                if local_name(attr.key.as_ref()) == b"activeObjective" {
                                    active_objective_id = Some(String::from_utf8_lossy(&attr.value).into_owned());
                                }
                            }
                        }
                        b"objective" => {
                            let mut oid = String::new();
                            for attr in e.attributes().flatten() {
                                if local_name(attr.key.as_ref()) == b"id" {
                                    oid = String::from_utf8_lossy(&attr.value).into_owned();
                                }
                            }
                            if !oid.is_empty() {
                                objectives.entry(oid.clone()).or_default();
                                current_objective_id = Some(oid);
                            }
                        }
                        b"fluxObjective" => {
                            let mut rxn_ref = String::new(); let mut coeff = 1.0;
                            for attr in e.attributes().flatten() {
                                let local = local_name(attr.key.as_ref());
                                let val = String::from_utf8_lossy(&attr.value);
                                match local {
                                    b"reaction" => rxn_ref = val.into_owned(),
                                    b"coefficient" => if let Ok(c) = val.parse::<f64>() { coeff = c; },
                                    _ => {}
                                }
                            }
                            if !rxn_ref.is_empty() {
                                if let Some(ref oid) = current_objective_id {
                                    if let Some(obj) = objectives.get_mut(oid) {
                                        obj.reaction_coefficients.push((rxn_ref, coeff));
                                    }
                                }
                            }
                        }
                        _ => {}
                    }
                }
                Ok(Event::End(e)) => {
                    if local_name(e.name().as_ref()) == b"objective" { current_objective_id = None; }
                }
                Ok(Event::Eof) => break,
                Err(e) => return Err(anyhow::anyhow!("XML pass 1: {}", e)),
                _ => {}
            }
            buf.clear();
        }
    }

    if let Some(bio_id) = determine_fbc_biomass(&active_objective_id, &objectives) {
        model.biomass_reaction = Some(bio_id);
    }

    // ── Pass 2: species, reactions, gene associations ──
    let mut current_species: Option<Metabolite> = None;
    let mut current_reaction: Option<Reaction> = None;
    let mut in_list_of_reactants = false;
    let mut in_list_of_products = false;
    let mut in_species_reference = false;
    let mut pending_species_ref: Option<(String, f64)> = None;
    let mut in_gene_assoc = false;
    let mut current_gene_ids: Vec<String> = Vec::new();
    let mut gene_assoc_stack: Vec<GeneAssocFrame> = Vec::new();
    let mut gene_assoc_root_children: Vec<String> = Vec::new();
    let mut gene_assoc_result: Option<String> = None;

    {
        let mut xr = Reader::from_str(&content);
        xr.trim_text(true);
        let mut buf = Vec::new();
        loop {
            match xr.read_event_into(&mut buf) {
                Ok(Event::Start(e)) => {
                    let tag_ref = e.name().as_ref().to_vec();
                    let tag_local = local_name(&tag_ref);
                    match tag_local {
                        b"species" => {
                            current_species = Some(parse_species_attrs(&e));
                        }
                        b"reaction" => {
                            let (rxn, lb_ref, ub_ref) = parse_reaction_attrs(&e, &parameters);
                            // Biomass fallback: name search
                            if model.biomass_reaction.is_none() {
                                let il = rxn.id.to_lowercase();
                                let nl = rxn.name.to_lowercase();
                                if nl.contains("biomass") || il.contains("biomass") || il.contains("growth") {
                                    model.biomass_reaction = Some(rxn.id.clone());
                                }
                            }
                            let _ = (lb_ref, ub_ref); // already resolved inside parse_reaction_attrs
                            current_reaction = Some(rxn);
                            current_gene_ids.clear();
                            gene_assoc_stack.clear();
                            gene_assoc_root_children.clear();
                            gene_assoc_result = None;
                        }
                        b"listOfReactants" => in_list_of_reactants = true,
                        b"listOfProducts" => in_list_of_products = true,
                        b"speciesReference" => {
                            let (sp, st) = parse_species_ref_attrs(&e);
                            in_species_reference = true;
                            pending_species_ref = Some((sp, st));
                        }
                        b"geneProductAssociation" => {
                            in_gene_assoc = true;
                            gene_assoc_stack.clear();
                            gene_assoc_root_children.clear();
                            gene_assoc_result = None;
                        }
                        b"and" if in_gene_assoc => gene_assoc_stack.push(GeneAssocFrame {
                            operator: "and".to_string(), children: Vec::new() }),
                        b"or" if in_gene_assoc => gene_assoc_stack.push(GeneAssocFrame {
                            operator: "or".to_string(), children: Vec::new() }),
                        _ => {}
                    }
                }
                Ok(Event::Empty(e)) => {
                    let tag_ref = e.name().as_ref().to_vec();
                    let tag_local = local_name(&tag_ref);
                    match tag_local {
                        b"species" => {
                            let m = parse_species_attrs(&e);
                            if !m.id.is_empty() { model.metabolites.insert(m.id.clone(), m); }
                        }
                        b"speciesReference" => {
                            let (sp, st) = parse_species_ref_attrs(&e);
                            add_species_ref_to_reaction(&mut current_reaction, &sp, st,
                                in_list_of_reactants, in_list_of_products);
                        }
                        b"geneProductRef" if in_gene_assoc => {
                            handle_gene_product_ref(&e, &model.genes,
                                &mut current_gene_ids, &mut gene_assoc_stack, &mut gene_assoc_root_children);
                        }
                        _ => {}
                    }
                }
                Ok(Event::End(e)) => {
                    let tag_ref = e.name().as_ref().to_vec();
                    let tag_local = local_name(&tag_ref);
                    match tag_local {
                        b"species" => {
                            if let Some(s) = current_species.take() {
                                model.metabolites.insert(s.id.clone(), s);
                            }
                        }
                        b"reaction" => {
                            if let Some(mut r) = current_reaction.take() {
                                current_gene_ids.sort(); current_gene_ids.dedup();
                                r.gene_ids = current_gene_ids.clone();
                                r.gene_reaction_rule = gene_assoc_result.take().unwrap_or_default();
                                model.reactions.insert(r.id.clone(), r);
                                current_gene_ids.clear();
                                gene_assoc_root_children.clear();
                            }
                        }
                        b"listOfReactants" => in_list_of_reactants = false,
                        b"listOfProducts" => in_list_of_products = false,
                        b"speciesReference" => {
                            if in_species_reference {
                                if let Some((ref sp, st)) = pending_species_ref {
                                    add_species_ref_to_reaction(&mut current_reaction, sp, st,
                                        in_list_of_reactants, in_list_of_products);
                                }
                                in_species_reference = false;
                                pending_species_ref = None;
                            }
                        }
                        b"geneProductAssociation" => {
                            in_gene_assoc = false;
                            if !gene_assoc_root_children.is_empty() {
                                gene_assoc_result = Some(if gene_assoc_root_children.len() == 1 {
                                    gene_assoc_root_children[0].clone()
                                } else { gene_assoc_root_children.join(" and ") });
                            }
                        }
                        b"and" | b"or" if in_gene_assoc => {
                            if let Some(frame) = gene_assoc_stack.pop() {
                                let sep = if frame.operator == "and" { " and " } else { " or " };
                                let expr = if frame.children.len() == 1 { frame.children[0].clone() }
                                    else if frame.children.is_empty() { String::new() }
                                    else { format!("({})", frame.children.join(sep)) };
                                if let Some(parent) = gene_assoc_stack.last_mut() {
                                    parent.children.push(expr);
                                } else { gene_assoc_root_children.push(expr); }
                            }
                        }
                        _ => {}
                    }
                }
                Ok(Event::Eof) => break,
                Err(e) => return Err(anyhow::anyhow!("XML pass 2: {}", e)),
                _ => {}
            }
            buf.clear();
        }
    }

    // ── Boundary EX reactions: O(n_rxns) indexed lookup, not O(n_b × n_rxns) ──
    let mut met_to_ex: HashMap<String, ()> = HashMap::new();
    for r in model.reactions.values() {
        if r.id.starts_with("R_EX_") || r.id.starts_with("EX_") {
            for (m, _) in &r.metabolites { met_to_ex.insert(m.clone(), ()); }
        }
    }
    let boundary_met_ids: Vec<String> = model.metabolites.values()
        .filter(|m| m.boundary).map(|m| m.id.clone()).collect();
    let mut added_boundary_ex = 0usize;
    for met_id in &boundary_met_ids {
        if !met_to_ex.contains_key(met_id) {
            let ex_id = format!("R_EX_boundary_{}", met_id);
            model.reactions.insert(ex_id.clone(), Reaction {
                id: ex_id, name: format!("Boundary exchange {}", met_id),
                metabolites: vec![(met_id.clone(), -1.0)],
                reversible: true, lower_bound: -1000.0, upper_bound: 1000.0,
                gene_reaction_rule: String::new(), gene_ids: Vec::new(),
            });
            added_boundary_ex += 1;
        }
    }

    let n_boundary = boundary_met_ids.len();
    let gene_linked = model.reactions.values().filter(|r| !r.gene_ids.is_empty()).count();
    let n_exchange = model.reactions.values()
        .filter(|r| r.id.starts_with("R_EX_") || r.id.starts_with("EX_")).count();
    eprintln!(
        "  SBML: {} parameters, {} genes, {} objectives, {} metabolites \
         ({} boundary, {} EX added), {} reactions ({} exchange, {} with genes), biomass: {:?}",
        parameters.len(), model.genes.len(), objectives.len(),
        model.metabolites.len(), n_boundary, added_boundary_ex,
        model.reactions.len(), n_exchange, gene_linked, model.biomass_reaction,
    );

    Ok(model)
}

// ── Helpers ──────────────────────────────────────────────────────────

fn local_name(tag: &[u8]) -> &[u8] {
    match tag.iter().position(|&b| b == b':') {
        Some(pos) => &tag[pos + 1..], None => tag,
    }
}

fn parse_species_attrs(e: &quick_xml::events::BytesStart) -> Metabolite {
    let (mut id, mut name, mut compartment, mut boundary, mut formula)
        = (String::new(), None, None, false, None);
    for attr in e.attributes().flatten() {
        // Match on the LOCAL name so namespace-prefixed attrs (e.g.
        // `fbc:chemicalFormula`, `fbc:charge`) are picked up regardless of
        // the document's namespace-prefix choice.
        match local_name(attr.key.as_ref()) {
            b"id"                => id = String::from_utf8_lossy(&attr.value).into_owned(),
            b"name"              => name = Some(String::from_utf8_lossy(&attr.value).into_owned()),
            b"compartment"       => compartment = Some(String::from_utf8_lossy(&attr.value).into_owned()),
            b"boundaryCondition" => boundary = &*attr.value == b"true",
            b"chemicalFormula"   => {
                let s = String::from_utf8_lossy(&attr.value).into_owned();
                if !s.trim().is_empty() { formula = Some(s); }
            }
            _ => {}
        }
    }
    Metabolite { id, name: name.unwrap_or_default(), compartment, boundary, formula }
}

fn parse_species_ref_attrs(e: &quick_xml::events::BytesStart) -> (String, f64) {
    let mut species_id = String::new(); let mut stoich = 1.0;
    for attr in e.attributes().flatten() {
        match attr.key.as_ref() {
            b"species" => species_id = String::from_utf8_lossy(&attr.value).into_owned(),
            b"stoichiometry" => {
                if let Ok(s) = std::str::from_utf8(&attr.value).unwrap_or("").parse::<f64>() {
                    stoich = s;
                }
            }
            _ => {}
        }
    }
    (species_id, stoich)
}

fn parse_reaction_attrs(
    e: &quick_xml::events::BytesStart, parameters: &HashMap<String, f64>,
) -> (Reaction, Option<String>, Option<String>) {
    let (mut id, mut name, mut reversible) = (String::new(), None, true);
    let (mut lb_ref, mut ub_ref) = (None, None);
    for attr in e.attributes().flatten() {
        let key = attr.key.as_ref();
        let local = local_name(key);
        let val = String::from_utf8_lossy(&attr.value);
        match local {
            b"id" => id = val.into_owned(),
            b"name" => name = Some(val.into_owned()),
            b"reversible" => reversible = val.as_ref() != "false",
            b"lowerFluxBound" => lb_ref = Some(val.into_owned()),
            b"upperFluxBound" => ub_ref = Some(val.into_owned()),
            _ => {}
        }
    }
    let lb = lb_ref.as_ref()
        .and_then(|r| parameters.get(r).copied().or_else(|| r.parse().ok()))
        .unwrap_or(if reversible { -1000.0 } else { 0.0 });
    let ub = ub_ref.as_ref()
        .and_then(|r| parameters.get(r).copied().or_else(|| r.parse().ok()))
        .unwrap_or(1000.0);
    let r = Reaction {
        id, name: name.unwrap_or_default(),
        metabolites: Vec::new(), reversible,
        lower_bound: lb, upper_bound: ub,
        gene_reaction_rule: String::new(), gene_ids: Vec::new(),
    };
    (r, lb_ref, ub_ref)
}

fn determine_fbc_biomass(
    active_id: &Option<String>, objectives: &HashMap<String, FbcObjective>,
) -> Option<String> {
    if let Some(ref aid) = active_id {
        if let Some(obj) = objectives.get(aid) {
            if obj.reaction_coefficients.len() == 1 {
                let (ref rxn_id, coeff) = obj.reaction_coefficients[0];
                if coeff > 0.0 { return Some(rxn_id.clone()); }
            }
            for (rxn_id, coeff) in &obj.reaction_coefficients {
                if *coeff > 0.0 { return Some(rxn_id.clone()); }
            }
        }
    }
    if objectives.len() == 1 {
        if let Some(obj) = objectives.values().next() {
            for (rxn_id, coeff) in &obj.reaction_coefficients {
                if *coeff > 0.0 { return Some(rxn_id.clone()); }
            }
        }
    }
    None
}

fn add_species_ref_to_reaction(
    current_reaction: &mut Option<Reaction>, species_id: &str, stoich: f64,
    in_reactants: bool, in_products: bool,
) {
    if species_id.is_empty() { return; }
    if let Some(ref mut reaction) = current_reaction {
        let coeff = if in_reactants { -stoich }
            else if in_products { stoich } else { stoich };
        reaction.metabolites.push((species_id.to_string(), coeff));
    }
}

fn handle_gene_product_ref(
    e: &quick_xml::events::BytesStart, genes: &HashMap<String, Gene>,
    current_gene_ids: &mut Vec<String>,
    gene_assoc_stack: &mut Vec<GeneAssocFrame>,
    gene_assoc_root_children: &mut Vec<String>,
) {
    for attr in e.attributes().flatten() {
        if local_name(attr.key.as_ref()) == b"geneProduct" {
            let gp = String::from_utf8_lossy(&attr.value).into_owned();
            if !gp.is_empty() {
                current_gene_ids.push(gp.clone());
                let label = genes.get(&gp).map(|g|
                    if g.label.is_empty() { g.id.clone() } else { g.label.clone() }
                ).unwrap_or_else(|| gp.clone());
                if let Some(frame) = gene_assoc_stack.last_mut() {
                    frame.children.push(label);
                } else {
                    gene_assoc_root_children.push(label);
                }
            }
        }
    }
}