//! Strict STL syntax and explicit facet-to-solid admission. STL carries no units:
//! callers must supply the source-length-to-inch factor. No sewing, normal
//! reversal, mesh repair, color guessing or smooth-surface reconstruction.
use crate::{
    geom::{Curve, Surface, Vec3},
    mass,
    neutral::{Body, Edge, Face},
    topology, KernelError,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Read,
    path::Path,
};
type Result<T> = std::result::Result<T, KernelError>;
fn invalid<T>(s: impl Into<String>) -> Result<T> {
    Err(KernelError::InvalidInput(s.into()))
}
fn unsupported<T>(s: impl Into<String>) -> Result<T> {
    Err(KernelError::Unsupported(s.into()))
}
const MAX_BYTES: usize = 64 * 1024 * 1024;
const MAX_FACETS: usize = 100_000;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Encoding {
    Ascii,
    Binary,
}
#[derive(Clone, Debug)]
pub struct Facet {
    pub normal: Vec3,
    pub vertices: [Vec3; 3],
}
#[derive(Clone, Debug)]
pub struct Mesh {
    pub encoding: Encoding,
    pub facets: Vec<Facet>,
}
fn facet(normal: Vec3, vertices: [Vec3; 3], scale: f64) -> Result<Facet> {
    let vertices = vertices.map(|p| p.scale(scale));
    if !normal.is_finite() || vertices.iter().any(|p| !p.is_finite()) {
        return invalid("STL contains a nonfinite normal/coordinate or unit conversion overflow");
    }
    let [a, b, c] = vertices;
    let n = b.sub(a).cross(c.sub(a));
    let length = n.length();
    if !length.is_finite() || length <= 1e-14 {
        return invalid("STL facet is degenerate or below the numerical area domain");
    }
    // An explicitly zero normal is unspecified. Winding remains authoritative;
    // nonzero supplied normals must agree rather than being silently replaced.
    if normal.length() != 0.
        && ((normal.length() - 1.).abs() > 1e-4
            || normal
                .scale(1. / normal.length())
                .distance(n.scale(1. / length))
                > 1e-4)
    {
        return invalid("STL facet normal disagrees with its vertex winding or is not unit length");
    }
    Ok(Facet { normal, vertices })
}
pub fn read_file(path: &Path, source_length_to_inches: f64) -> Result<Mesh> {
    let file = std::fs::File::open(path)
        .map_err(|e| KernelError::InvalidInput(format!("cannot read STL: {e}")))?;
    let mut bytes = vec![];
    file.take(MAX_BYTES as u64 + 1)
        .read_to_end(&mut bytes)
        .map_err(|e| KernelError::InvalidInput(format!("cannot read STL: {e}")))?;
    read(&bytes, source_length_to_inches)
}
pub fn read(bytes: &[u8], source_length_to_inches: f64) -> Result<Mesh> {
    let scale = source_length_to_inches;
    if !scale.is_finite() || scale <= 0. {
        return invalid("STL requires an explicit finite positive source-length-to-inch factor");
    }
    if bytes.len() > MAX_BYTES {
        return unsupported("STL exceeds the 64 MiB input budget");
    }
    // Exact binary record length identifies binary files even if their arbitrary
    // 80-byte header begins with 'solid'. Classification precedes parsing.
    let binary_count = bytes
        .get(80..84)
        .map(|v| u32::from_le_bytes(v.try_into().unwrap()) as usize);
    if binary_count.is_some_and(|n| 84u64 + 50 * n as u64 == bytes.len() as u64) {
        let n = binary_count.unwrap();
        if n == 0 || n > MAX_FACETS {
            return unsupported("STL requires 1..100000 facets");
        }
        let mut facets = Vec::with_capacity(n);
        for (i, chunk) in bytes[84..].chunks_exact(50).enumerate() {
            if chunk[48] != 0 || chunk[49] != 0 {
                return unsupported(format!(
                    "STL facet {i} has nonzero attributes; color conventions are not inferred"
                ));
            }
            let number = |at| f32::from_le_bytes(chunk[at..at + 4].try_into().unwrap()) as f64;
            let point = |at| Vec3::new(number(at), number(at + 4), number(at + 8));
            facets.push(facet(point(0), [point(12), point(24), point(36)], scale)?);
        }
        return Ok(Mesh {
            encoding: Encoding::Binary,
            facets,
        });
    }
    if !bytes.is_ascii() || bytes.contains(&0) {
        return invalid("STL binary count/length mismatch or non-ASCII text");
    }
    let text = std::str::from_utf8(bytes).unwrap();
    let mut lines = text.lines().map(str::trim).filter(|s| !s.is_empty());
    let label = |line: &str, key: &str| {
        line.strip_prefix(key)
            .filter(|tail| tail.is_empty() || tail.starts_with(char::is_whitespace))
            .map(|s| s.trim().to_owned())
    };
    let Some(name) = lines.next().and_then(|line| label(line, "solid")) else {
        return invalid(
            "STL text must begin with a solid record (or binary count/length is invalid)",
        );
    };
    let point = |line: &str, prefix: &str| -> Result<Vec3> {
        let Some(tail) = line.strip_prefix(prefix) else {
            return invalid(format!("STL expected {prefix}"));
        };
        if !tail.starts_with(char::is_whitespace) {
            return invalid("malformed STL vector record");
        }
        let mut tokens = tail.split_whitespace();
        let mut values = [0.; 3];
        for value in &mut values {
            let token = tokens.next().ok_or_else(|| {
                KernelError::InvalidInput("STL vector needs three coordinates".into())
            })?;
            *value = token
                .parse::<f64>()
                .map_err(|_| KernelError::InvalidInput("invalid STL number".into()))?;
        }
        if tokens.next().is_some() {
            return invalid("STL vector requires exactly three coordinates");
        }
        Ok(Vec3::new(values[0], values[1], values[2]))
    };
    let mut facets = vec![];
    while let Some(line) = lines.next() {
        if let Some(end) = label(line, "endsolid") {
            if !end.is_empty() && end != name {
                return invalid("STL endsolid name differs from solid name");
            }
            if lines.next().is_some() {
                return unsupported(
                    "STL trailing records or multiple ASCII solids are not supported",
                );
            }
            if facets.is_empty() {
                return invalid("STL has no facets");
            }
            return Ok(Mesh {
                encoding: Encoding::Ascii,
                facets,
            });
        }
        if facets.len() >= MAX_FACETS {
            return unsupported("STL exceeds 100000 facets");
        }
        let normal = point(line, "facet normal")?;
        if lines.next() != Some("outer loop") {
            return invalid("STL expected outer loop");
        }
        let mut vertices = [Vec3::new(0., 0., 0.); 3];
        for v in &mut vertices {
            *v = point(
                lines
                    .next()
                    .ok_or_else(|| KernelError::InvalidInput("truncated STL vertices".into()))?,
                "vertex",
            )?;
        }
        if lines.next() != Some("endloop") || lines.next() != Some("endfacet") {
            return invalid("STL expected endloop/endfacet after exactly three vertices");
        }
        facets.push(facet(normal, vertices, scale)?);
    }
    invalid("STL is missing endsolid")
}

/// Separating axes certify disjoint interiors. When triangles share topology,
/// contact is permitted only if the supporting plane confines at least one
/// triangle's entire contact set to that exact shared vertex/edge.
fn separated(a: [Vec3; 3], b: [Vec3; 3], ia: [usize; 3], ib: [usize; 3]) -> bool {
    let origin = a[0];
    let a = a.map(|p| p.sub(origin));
    let b = b.map(|p| p.sub(origin));
    let ea = std::array::from_fn::<_, 3, _>(|i| a[(i + 1) % 3].sub(a[i]));
    let eb = std::array::from_fn::<_, 3, _>(|i| b[(i + 1) % 3].sub(b[i]));
    let (na, nb) = (ea[0].cross(ea[1]), eb[0].cross(eb[1]));
    let mut axes = vec![na, nb];
    for e in ea {
        axes.push(e.cross(na));
        for f in eb {
            axes.push(e.cross(f));
        }
    }
    for e in eb {
        axes.push(e.cross(nb));
    }
    let round = 512. * f64::EPSILON * a.iter().chain(&b).map(|p| p.length()).fold(1., f64::max);
    for axis in axes {
        let len = axis.length();
        if len == 0. || !len.is_finite() {
            continue;
        }
        let axis = axis.scale(1. / len);
        let pa = a.map(|p| p.dot(axis));
        let pb = b.map(|p| p.dot(axis));
        let lo = |p: [f64; 3]| p.into_iter().fold(f64::INFINITY, f64::min);
        let hi = |p: [f64; 3]| p.into_iter().fold(f64::NEG_INFINITY, f64::max);
        for (x, y, ix, iy) in [(pa, pb, ia, ib), (pb, pa, ib, ia)] {
            let (h, l) = (hi(x), lo(y));
            if h + round < l {
                return true;
            }
            if (h - l).abs() <= round && ix.iter().any(|i| iy.contains(i)) {
                let only_shared_x = (0..3).all(|k| iy.contains(&ix[k]) || x[k] < h - round);
                let only_shared_y = (0..3).all(|k| ix.contains(&iy[k]) || y[k] > l + round);
                if only_shared_x || only_shared_y {
                    return true;
                }
            }
        }
    }
    false
}

pub fn to_solid(mesh: &Mesh, tolerance: f64) -> Result<Body> {
    if !tolerance.is_finite()
        || tolerance <= 0.
        || mesh.facets.is_empty()
        || mesh.facets.len() > MAX_FACETS
    {
        return invalid("invalid STL solid tolerance/facet count");
    }
    let mut body = Body::default();
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    let mut triangles = vec![];
    let mut owners: Vec<Vec<usize>> = vec![];
    for (f, source) in mesh.facets.iter().enumerate() {
        let source = facet(source.normal, source.vertices, 1.)?;
        let mut ids = [0; 3];
        for i in 0..3 {
            let p = source.vertices[i];
            let bits = |x: f64| if x == 0. { 0 } else { x.to_bits() };
            ids[i] = *vertices
                .entry([bits(p.x), bits(p.y), bits(p.z)])
                .or_insert_with(|| {
                    let i = body.vertices.len();
                    body.vertices.push(p);
                    i
                });
        }
        let mut lp = vec![];
        for i in 0..3 {
            let (a, b) = (ids[i], ids[(i + 1) % 3]);
            let key = (a.min(b), a.max(b));
            let e = *edges.entry(key).or_insert_with(|| {
                let e = body.edges.len();
                body.edges.push(Edge {
                    curve: Curve::Line {
                        a: body.vertices[key.0],
                        b: body.vertices[key.1],
                    },
                    start: key.0,
                    end: key.1,
                });
                owners.push(vec![]);
                e
            });
            owners[e].push(f);
            lp.push((e, a == key.0));
        }
        let [a, b, c] = source.vertices;
        body.faces.push(Face {
            surface: Surface::Plane {
                p: a,
                n: b.sub(a).cross(c.sub(a)).unit()?,
                x: b.sub(a).unit()?,
            },
            sense: true,
            loops: vec![lp],
        });
        triangles.push(ids);
    }
    topology::check(&body, tolerance)?;
    // Every vertex link must be one cycle; edge counts alone miss pinched shells.
    let mut links: Vec<BTreeMap<usize, Vec<usize>>> = vec![BTreeMap::new(); body.vertices.len()];
    for ids in &triangles {
        for i in 0..3 {
            let (v, a, b) = (ids[i], ids[(i + 1) % 3], ids[(i + 2) % 3]);
            links[v].entry(a).or_default().push(b);
            links[v].entry(b).or_default().push(a);
        }
    }
    for link in links {
        if link.values().any(|n| n.len() != 2) {
            return invalid("STL solid has a nonmanifold vertex link");
        }
        let mut seen = BTreeSet::new();
        let mut todo = vec![*link.keys().next().unwrap()];
        while let Some(v) = todo.pop() {
            if seen.insert(v) {
                todo.extend(&link[&v]);
            }
        }
        if seen.len() != link.len() {
            return invalid("STL solid has a pinched vertex");
        }
    }
    let mut seen = BTreeSet::new();
    let mut todo = vec![0];
    while let Some(f) = todo.pop() {
        if seen.insert(f) {
            for &(e, _) in &body.faces[f].loops[0] {
                todo.extend(&owners[e]);
            }
        }
    }
    if seen.len() != body.faces.len() {
        return unsupported("STL solid conversion requires one connected shell; use mesh intake for disconnected facets");
    }
    // Sweep bounding boxes before the bounded exact-feature intersection audit.
    let bounds: Vec<_> = mesh
        .facets
        .iter()
        .map(|f| {
            let p = f.vertices.map(|p| [p.x, p.y, p.z]);
            (
                std::array::from_fn::<_, 3, _>(|k| {
                    p.iter().map(|p| p[k]).fold(f64::INFINITY, f64::min)
                }),
                std::array::from_fn::<_, 3, _>(|k| {
                    p.iter().map(|p| p[k]).fold(f64::NEG_INFINITY, f64::max)
                }),
            )
        })
        .collect();
    let mut order: Vec<_> = (0..triangles.len()).collect();
    order.sort_by(|&i, &j| bounds[i].0[0].total_cmp(&bounds[j].0[0]));
    let mut checks = 0;
    for x in 0..order.len() {
        let i = order[x];
        for &j in &order[x + 1..] {
            if bounds[j].0[0] > bounds[i].1[0] {
                break;
            }
            checks += 1;
            if checks > 2_000_000 {
                return unsupported("STL intersection audit exceeds 2000000 candidate pairs");
            }
            if (1..3).any(|k| bounds[i].1[k] < bounds[j].0[k] || bounds[j].1[k] < bounds[i].0[k]) {
                continue;
            }
            if !separated(
                mesh.facets[i].vertices,
                mesh.facets[j].vertices,
                triangles[i],
                triangles[j],
            ) {
                return invalid(format!(
                    "STL facets {i} and {j} intersect, overlap, or have unresolved contact"
                ));
            }
        }
    }
    let properties = mass::measure(&body, tolerance)?;
    if properties.volume <= 0. {
        return invalid("STL shell is inward oriented or has zero volume");
    }
    Ok(body)
}

/// Op 18: [1, mode (0 mesh / 1 closed faceted solid), explicit unit factor,
/// UTF8 path length, bytes]. No repair and no source/destination writes.
pub fn read_op(data: &[f64], tolerance: f64) -> Result<Vec<f64>> {
    let mut r = crate::neutral::Reader::with_tolerance(data, tolerance);
    if r.next()? != 1. {
        return invalid("invalid STL envelope version");
    }
    let solid = r.flag()?;
    let scale = r.next()?;
    let n = r.count()?;
    if n == 0 || n > 32768 {
        return invalid("invalid STL path length");
    }
    let bytes = (0..n)
        .map(|_| r.index(256, "STL path byte").map(|v| v as u8))
        .collect::<Result<Vec<_>>>()?;
    r.finish()?;
    if bytes.contains(&0) {
        return invalid("STL path contains a zero byte");
    }
    let path = String::from_utf8(bytes)
        .map_err(|_| KernelError::InvalidInput("STL path is not UTF-8".into()))?;
    let mesh = read_file(Path::new(&path), scale)?;
    if solid {
        let body = to_solid(&mesh, tolerance)?;
        let volume = mass::measure(&body, tolerance)?.volume;
        let data = crate::neutral::write_body(&body);
        let mut out = vec![1., 1., volume, data.len() as f64];
        out.extend(data);
        Ok(out)
    } else {
        let mut out = vec![1., 0., mesh.facets.len() as f64];
        for f in mesh.facets {
            for p in [f.normal, f.vertices[0], f.vertices[1], f.vertices[2]] {
                out.extend([p.x, p.y, p.z]);
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn separation_rejects_interior_crossings_and_accepts_only_shared_contact() {
        let p = |x, y, z| Vec3::new(x, y, z);
        let a = [p(0., 0., 0.), p(2., 0., 0.), p(0., 2., 0.)];
        let coplanar = [p(0.2, 0.2, 0.), p(1., 0.2, 0.), p(0.2, 1., 0.)];
        let crossing = [p(0.2, 0.2, -1.), p(0.2, 0.2, 1.), p(1., 0.2, 0.)];
        for b in [coplanar, crossing] {
            assert!(!separated(a, b, [0, 1, 2], [3, 4, 5]));
        }
        let adjacent = [a[1], a[0], p(1., -1., 0.)];
        assert!(separated(a, adjacent, [0, 1, 2], [1, 0, 3]));
        assert!(!separated(a, adjacent, [0, 1, 2], [3, 4, 5]));
        let overlap = [a[1], a[0], p(0.5, 0.5, 0.)];
        assert!(!separated(a, overlap, [0, 1, 2], [1, 0, 3]));
    }
}
