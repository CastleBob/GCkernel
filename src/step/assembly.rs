//! Explicit assembly occurrences. Definition bodies are cached in inches;
//! each occurrence carries its source path and a proper rigid world placement.
use super::{context, document, invalid, syntax, unsupported, ImportedBody, Result};
use crate::{
    frame::{Frame, Point3, Rigid},
    mass, neutral, topology,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Read,
    path::Path,
};
#[derive(Clone, Debug)]
pub struct Instance {
    /// Root representation followed by usage/identity-relationship IDs.
    pub path: Vec<u64>,
    pub representation: u64,
    pub placement: Rigid,
    /// Detached geometry already placed in world coordinates, in inches.
    pub item: ImportedBody,
}
struct Rep {
    scale: f64,
    angle: f64,
    items: Vec<u64>,
    solids: Vec<u64>,
}
struct Link {
    child: u64,
    id: u64,
    motion: Rigid,
}
pub fn read_file(path: &Path, tolerance: f64) -> Result<Vec<Instance>> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)
        .map_err(|e| crate::KernelError::InvalidInput(e.to_string()))?
        .take(syntax::MAX_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| crate::KernelError::InvalidInput(e.to_string()))?;
    read(&bytes, tolerance)
}
pub fn read(bytes: &[u8], tolerance: f64) -> Result<Vec<Instance>> {
    let doc = document(bytes, tolerance)?;
    let mut reps = BTreeMap::new();
    for (&id, e) in &doc.exchange.entities {
        for kind in e.components.keys() {
            if kind.ends_with("SHAPE_REPRESENTATION")
                && !matches!(
                    kind.as_str(),
                    "SHAPE_REPRESENTATION"
                        | "ADVANCED_BREP_SHAPE_REPRESENTATION"
                        | "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"
                )
            {
                return unsupported(format!(
                    "STEP #{id} unsupported assembly representation {kind}"
                ));
            }
        }
        for kind in [
            "MAPPED_ITEM",
            "BREP_WITH_VOIDS",
            "FACETED_BREP",
            "SHELL_BASED_SURFACE_MODEL",
            "CSG_SOLID",
            "SOLID_REPLICA",
            "EXTERNALLY_DEFINED_REPRESENTATION",
        ] {
            if e.has(kind) {
                return unsupported(format!("STEP #{id} unsupported assembly item {kind}"));
            }
        }
        for kind in ["SHAPE_REPRESENTATION", "ADVANCED_BREP_SHAPE_REPRESENTATION"] {
            if !e.has(kind) {
                continue;
            }
            let a = e.args(kind, 3)?;
            let (scale, angle) = doc.units(a[2].reference()?)?;
            let items = a[1]
                .list()?
                .iter()
                .map(|v| v.reference())
                .collect::<Result<Vec<_>>>()?;
            let mut solids = Vec::new();
            let mut seen = BTreeSet::new();
            for &item in &items {
                if !seen.insert(item) {
                    return invalid("duplicate assembly representation item");
                }
                let e = doc.entity(item)?;
                if e.has("MANIFOLD_SOLID_BREP") {
                    solids.push(item);
                } else if e.has("AXIS2_PLACEMENT_3D") {
                    doc.frame(item)?;
                } else {
                    return unsupported(format!(
                        "STEP #{item} unsupported assembly representation item {}",
                        e.kind()
                    ));
                }
            }
            reps.insert(
                id,
                Rep {
                    scale,
                    angle,
                    items,
                    solids,
                },
            );
        }
    }
    let mut product_rep = BTreeMap::new();
    for e in doc
        .exchange
        .entities
        .values()
        .filter(|e| e.has("SHAPE_DEFINITION_REPRESENTATION"))
    {
        let a = e.args("SHAPE_DEFINITION_REPRESENTATION", 2)?;
        let product = doc.args(a[0].reference()?, "PRODUCT_DEFINITION_SHAPE", 3)?[2].reference()?;
        let rep = a[1].reference()?;
        if !reps.contains_key(&rep) {
            return unsupported("product shape does not identify a supported representation");
        }
        if product_rep.insert(product, rep).is_some() {
            return unsupported("product has multiple shape definitions");
        }
    }
    let mut links: BTreeMap<u64, Vec<Link>> = BTreeMap::new();
    let mut used_relationships = BTreeSet::new();
    let mut used_usages = BTreeSet::new();
    for e in doc
        .exchange
        .entities
        .values()
        .filter(|e| e.has("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION"))
    {
        let a = e.args("CONTEXT_DEPENDENT_SHAPE_REPRESENTATION", 2)?;
        let relation = a[0].reference()?;
        let usage = doc.args(a[1].reference()?, "PRODUCT_DEFINITION_SHAPE", 3)?[2].reference()?;
        if !used_usages.insert(usage) || !used_relationships.insert(relation) {
            return invalid("duplicate assembly usage/relationship");
        }
        let usage_args = doc.args(usage, "NEXT_ASSEMBLY_USAGE_OCCURRENCE", 6)?;
        let lookup = |id| {
            product_rep.get(&id).copied().ok_or_else(|| {
                crate::KernelError::InvalidInput("assembly product has no shape definition".into())
            })
        };
        let parent = lookup(usage_args[3].reference()?)?;
        let child = lookup(usage_args[4].reference()?)?;
        let r = doc.entity(relation)?;
        let (a, transform) = if r.has("REPRESENTATION_RELATIONSHIP") {
            (
                r.args("REPRESENTATION_RELATIONSHIP", 4)?,
                r.args("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION", 1)?[0].reference()?,
            )
        } else {
            let a = r.args("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION", 5)?;
            (&a[..4], a[4].reference()?)
        };
        let (one, two) = (a[2].reference()?, a[3].reference()?);
        if !((one == parent && two == child) || (one == child && two == parent)) || parent == child
        {
            return invalid("assembly usage and representation relationship disagree");
        }
        let a = doc.args(transform, "ITEM_DEFINED_TRANSFORMATION", 4)?;
        let frame = |rep: u64, id: u64| -> Result<Frame> {
            if !reps[&rep].items.contains(&id) {
                return invalid("assembly placement is not an item of its representation");
            }
            let (p, n, x) = doc.frame(id)?;
            Frame::new(Point3::from_vec3(p.scale(reps[&rep].scale)), x, n.cross(x))
        };
        let motion = Rigid::from_frames(
            &frame(one, a[2].reference()?)?,
            &frame(two, a[3].reference()?)?,
        );
        let motion = if one == child {
            motion
        } else {
            motion.inverse()
        };
        links.entry(parent).or_default().push(Link {
            child,
            id: usage,
            motion,
        });
    }
    for (&id, e) in &doc.exchange.entities {
        if e.has("NEXT_ASSEMBLY_USAGE_OCCURRENCE") && !used_usages.contains(&id) {
            return invalid("assembly usage has no context-dependent placement");
        }
        if e.has("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
            && !used_relationships.contains(&id)
        {
            return unsupported("unbound representation transformation");
        }
        if e.has("SHAPE_REPRESENTATION_RELATIONSHIP")
            && !e.has("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION")
        {
            let a = e.args("SHAPE_REPRESENTATION_RELATIONSHIP", 4)?;
            let (parent, child) = (a[2].reference()?, a[3].reference()?);
            if !reps.contains_key(&parent) || !reps.contains_key(&child) || parent == child {
                return invalid("invalid shape representation relationship");
            }
            // Explicit identity relationship: units are already converted per
            // representation. No inferred placement or missing transformation.
            links.entry(parent).or_default().push(Link {
                child,
                id,
                motion: Rigid::IDENTITY,
            });
        }
    }
    // Validate the entire graph, including cycles disconnected from roots.
    fn acyclic(
        id: u64,
        links: &BTreeMap<u64, Vec<Link>>,
        active: &mut BTreeSet<u64>,
        done: &mut BTreeSet<u64>,
    ) -> Result<()> {
        if done.contains(&id) {
            return Ok(());
        }
        if active.len() >= 64 || !active.insert(id) {
            return invalid("cyclic or excessively deep assembly graph");
        }
        if let Some(children) = links.get(&id) {
            for link in children {
                acyclic(link.child, links, active, done)?;
            }
        }
        active.remove(&id);
        done.insert(id);
        Ok(())
    }
    let mut done = BTreeSet::new();
    for &id in reps.keys() {
        acyclic(id, &links, &mut BTreeSet::new(), &mut done)?;
    }
    let mut definitions = BTreeMap::new();
    for (&rep, r) in &reps {
        for &id in &r.solids {
            if definitions
                .insert(
                    (rep, id),
                    doc.body(id, r.scale, r.angle).map_err(|e| context(e, id))?,
                )
                .is_some()
            {
                return invalid("duplicate solid definition");
            }
        }
    }
    let defined: BTreeSet<_> = definitions.keys().map(|(_, id)| *id).collect();
    if defined.is_empty() {
        return invalid("assembly contains no closed solids");
    }
    for (&id, e) in &doc.exchange.entities {
        if e.has("MANIFOLD_SOLID_BREP") && !defined.contains(&id) {
            return invalid("assembly contains an unrepresented solid");
        }
    }
    let children: BTreeSet<_> = links.values().flatten().map(|l| l.child).collect();
    let mut pending = Vec::new();
    for &id in reps.keys().filter(|id| !children.contains(id)) {
        pending.push((id, vec![id], Rigid::IDENTITY));
    }
    let mut out = Vec::new();
    let mut visited = BTreeSet::new();
    let mut visits = 0;
    while let Some((id, path, motion)) = pending.pop() {
        visits += 1;
        if visits > 10000 || path.len() > 64 {
            return unsupported("assembly occurrence expansion exceeds budget");
        }
        for &solid in &reps[&id].solids {
            if out.len() >= 1024 {
                return unsupported("assembly exceeds 1024 solid occurrences");
            }
            let mut item = definitions[&(id, solid)].clone();
            item.body = motion.body(&item.body);
            let data = neutral::write_body(&item.body);
            let mut r = neutral::Reader::with_tolerance(&data, tolerance);
            item.body = r.body()?;
            r.finish()?;
            topology::check(&item.body, tolerance)?;
            item.volume = mass::measure(&item.body, tolerance)?.volume;
            visited.insert((id, solid));
            out.push(Instance {
                path: path.clone(),
                representation: id,
                placement: motion,
                item,
            });
        }
        if let Some(children) = links.get(&id) {
            for child in children.iter().rev() {
                let mut path = path.clone();
                path.push(child.id);
                pending.push((child.child, path, child.motion.then(&motion)));
            }
        }
    }
    if visited.len() != definitions.len() {
        return invalid("assembly contains an unreachable definition");
    }
    Ok(out)
}
