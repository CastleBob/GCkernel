//! Independently authored rectilinear solid subtraction. The input's actual
//! coordinates partition space into cells on which membership is constant.
//! This is an exact planar arrangement for the admitted domain, not voxel
//! sampling: no resolution, snapping, offsets, meshes or backup engine.
//! Every operand boundary is verified against the classified material. Results
//! retain analytic planar surfaces, oriented loops and shared edges/vertices.
use crate::geom::{Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError as CutError;
use std::collections::{BTreeMap, BTreeSet};

type Result<T> = std::result::Result<T, CutError>;
const MAX_COORDINATES: usize = 96;
const MAX_CELLS: usize = 200_000;
const MAX_CHECKS: usize = 20_000_000;
fn invalid(m: &str) -> CutError {
    CutError::InvalidInput(m.into())
}
fn unsupported(m: &str) -> CutError {
    CutError::Unsupported(m.into())
}
fn xyz(p: Vec3) -> [f64; 3] {
    [p.x, p.y, p.z]
}
fn vec(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}
fn finite(p: Vec3) -> bool {
    xyz(p).iter().all(|v| v.is_finite())
}

pub(super) fn cut(a: &Body, b: &Body, tolerance: f64) -> Result<Body> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(invalid("cut tolerance must be finite and positive"));
    }
    let af = validate_faces(a)?;
    let bf = validate_faces(b)?;
    if af.is_empty() && bf.is_empty() {
        return Ok(Body::default());
    }
    let grid = Grid::new(a, b, tolerance)?;
    let ac = grid.membership(&af)?;
    let bc = grid.membership(&bf)?;
    let cells: Vec<bool> = ac.iter().zip(&bc).map(|(a, b)| *a && !*b).collect();
    grid.require_connected_material(&cells)?;
    let result = grid.make_body(&cells);
    let faces = validate_faces(&result).map_err(|_| {
        unsupported("cut produces a nonmanifold boundary outside the supported result domain")
    })?;
    if grid.membership(&faces)? != cells {
        return Err(invalid(
            "result boundary disagrees with calculated material",
        ));
    }
    Ok(result)
}

#[derive(Debug)]
struct FlatFace {
    axis: usize,
    coordinate: f64,
    sign: i8,
    loops: Vec<Vec<[f64; 2]>>,
}

fn axis_of(d: Vec3) -> Result<(usize, i8)> {
    if !finite(d) {
        return Err(invalid("nonfinite direction"));
    }
    let a = xyz(d);
    let axes: Vec<usize> = (0..3).filter(|i| a[*i] != 0.).collect();
    if axes.len() != 1 {
        return Err(unsupported(
            "cut requires axis-aligned planar faces and straight edges",
        ));
    }
    Ok((axes[0], if a[axes[0]] > 0. { 1 } else { -1 }))
}

fn validate_faces(body: &Body) -> Result<Vec<FlatFace>> {
    if body.vertices.is_empty() && body.edges.is_empty() && body.faces.is_empty() {
        return Ok(Vec::new());
    }
    if body.faces.is_empty() || body.edges.is_empty() || body.vertices.is_empty() {
        return Err(invalid("incomplete body topology"));
    }
    if body.faces.len() > 20_000 || body.edges.len() > 80_000 || body.vertices.len() > 80_000 {
        return Err(unsupported("rectilinear topology limit exceeded"));
    }
    if body.vertices.iter().any(|p| !finite(*p)) {
        return Err(invalid("nonfinite vertex"));
    }
    for e in &body.edges {
        let a = *body
            .vertices
            .get(e.start)
            .ok_or_else(|| invalid("invalid edge start"))?;
        let b = *body
            .vertices
            .get(e.end)
            .ok_or_else(|| invalid("invalid edge end"))?;
        match &e.curve {
            Curve::Line { a: ca, b: cb } => {
                if xyz(*ca) != xyz(a) || xyz(*cb) != xyz(b) {
                    return Err(invalid("curve endpoints disagree with topology"));
                }
                if e.start == e.end || xyz(a) == xyz(b) {
                    return Err(invalid("zero-length edge"));
                }
                axis_of(b.sub(a))?;
            }
            _ => {
                return Err(unsupported(
                    "curved edges are not supported by rectilinear cut",
                ))
            }
        }
    }
    let mut uses = vec![[0usize; 2]; body.edges.len()];
    let mut edge_faces: Vec<Vec<usize>> = vec![Vec::new(); body.edges.len()];
    let mut faces = Vec::new();
    let mut checks = 0usize;
    for (face_index, f) in body.faces.iter().enumerate() {
        let (p, n, x) = match &f.surface {
            Surface::Plane { p, n, x } => (*p, *n, *x),
            _ => {
                return Err(unsupported(
                    "curved faces are not supported by rectilinear cut",
                ))
            }
        };
        if !finite(p) || !finite(x) {
            return Err(invalid("nonfinite plane frame"));
        }
        let (axis, normal_sign) = axis_of(n)?;
        if xyz(x)[axis] != 0. || x.dot(x) == 0. || !x.dot(x).is_finite() {
            return Err(invalid("invalid plane reference direction"));
        }
        let sign = normal_sign * if f.sense { 1 } else { -1 };
        let coordinate = xyz(p)[axis];
        if f.loops.is_empty() {
            return Err(invalid("face has no boundaries"));
        }
        let mut loops = Vec::new();
        for lp in &f.loops {
            if lp.len() < 4 {
                return Err(invalid("rectilinear loop needs four or more edges"));
            }
            checks = checks.saturating_add(lp.len().saturating_mul(lp.len()));
            if checks > MAX_CHECKS {
                return Err(unsupported("loop validation budget exceeded"));
            }
            let mut polygon = Vec::new();
            for (i, (index, forward)) in lp.iter().enumerate() {
                let e = body
                    .edges
                    .get(*index)
                    .ok_or_else(|| invalid("invalid loop edge"))?;
                uses[*index][usize::from(*forward)] += 1;
                edge_faces[*index].push(face_index);
                let (start, end) = if *forward {
                    (e.start, e.end)
                } else {
                    (e.end, e.start)
                };
                let (next_index, next_forward) = lp[(i + 1) % lp.len()];
                let next = body
                    .edges
                    .get(next_index)
                    .ok_or_else(|| invalid("invalid loop edge"))?;
                if end != if next_forward { next.start } else { next.end } {
                    return Err(invalid("disconnected face loop"));
                }
                let point = xyz(body.vertices[start]);
                if point[axis] != coordinate {
                    return Err(invalid("loop vertex is off plane"));
                }
                polygon.push([point[(axis + 1) % 3], point[(axis + 2) % 3]]);
            }
            validate_simple_loop(&polygon)?;
            loops.push(polygon);
        }
        for i in 0..loops.len() {
            for j in i + 1..loops.len() {
                checks = checks.saturating_add(loops[i].len().saturating_mul(loops[j].len()));
                if checks > MAX_CHECKS {
                    return Err(unsupported("face validation budget exceeded"));
                }
                if loops_intersect(&loops[i], &loops[j]) {
                    return Err(invalid("face loops cross or touch"));
                }
            }
            let depth = (0..loops.len())
                .filter(|j| *j != i && in_loop(loops[i][0], &loops[*j]))
                .count();
            let expected = sign as f64 * if depth % 2 == 0 { 1. } else { -1. };
            let area = signed_area(&loops[i]);
            if !area.is_finite() || area * expected <= 0. {
                return Err(invalid("loop orientation disagrees with face sense"));
            }
        }
        faces.push(FlatFace {
            axis,
            coordinate,
            sign,
            loops,
        });
    }
    if uses.iter().any(|u| *u != [1, 1]) {
        return Err(invalid(
            "every edge needs two opposite oriented occurrences",
        ));
    }
    // Each vertex must have a single connected fan of incident faces. Merely
    // checking edge occurrences misses shells touching at a lone vertex.
    let mut fans: Vec<Vec<(usize, usize)>> = vec![Vec::new(); body.vertices.len()];
    for (e, fs) in body.edges.iter().zip(&edge_faces) {
        fans[e.start].push((fs[0], fs[1]));
        fans[e.end].push((fs[0], fs[1]));
    }
    for fan in fans {
        if fan.is_empty() {
            return Err(invalid("unused vertex"));
        }
        let mut seen = BTreeSet::from([fan[0].0]);
        loop {
            let old = seen.len();
            for (a, b) in &fan {
                if seen.contains(a) || seen.contains(b) {
                    seen.insert(*a);
                    seen.insert(*b);
                }
            }
            if old == seen.len() {
                break;
            }
        }
        if fan
            .iter()
            .any(|(a, b)| !seen.contains(a) || !seen.contains(b))
        {
            return Err(invalid("nonmanifold vertex fan"));
        }
    }
    Ok(faces)
}

fn signed_area(p: &[[f64; 2]]) -> f64 {
    let o = p[0];
    (0..p.len())
        .map(|i| {
            let a = p[i];
            let b = p[(i + 1) % p.len()];
            (a[0] - o[0]) * (b[1] - o[1]) - (b[0] - o[0]) * (a[1] - o[1])
        })
        .sum::<f64>()
        * 0.5
}
fn segments_touch(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    (0..2).all(|i| a[i].min(b[i]).max(c[i].min(d[i])) <= a[i].max(b[i]).min(c[i].max(d[i])))
}
fn validate_simple_loop(p: &[[f64; 2]]) -> Result<()> {
    for i in 0..p.len() {
        let a = p[i];
        let b = p[(i + 1) % p.len()];
        let c = p[(i + 2) % p.len()];
        if (a[0] == b[0]) == (a[1] == b[1]) {
            return Err(invalid("invalid rectilinear loop segment"));
        }
        if (b[0] - a[0]) * (c[0] - b[0]) + (b[1] - a[1]) * (c[1] - b[1]) < 0. {
            return Err(invalid("loop doubles back"));
        }
        for j in i + 1..p.len() {
            if j == i + 1 || (i == 0 && j == p.len() - 1) {
                continue;
            }
            if segments_touch(a, b, p[j], p[(j + 1) % p.len()]) {
                return Err(invalid("self-intersecting loop"));
            }
        }
    }
    Ok(())
}
fn loops_intersect(a: &[[f64; 2]], b: &[[f64; 2]]) -> bool {
    (0..a.len()).any(|i| {
        (0..b.len()).any(|j| segments_touch(a[i], a[(i + 1) % a.len()], b[j], b[(j + 1) % b.len()]))
    })
}
fn in_loop(p: [f64; 2], lp: &[[f64; 2]]) -> bool {
    let mut inside = false;
    for i in 0..lp.len() {
        let a = lp[i];
        let b = lp[(i + 1) % lp.len()];
        if (a[1] > p[1]) != (b[1] > p[1]) && a[0] > p[0] {
            inside = !inside;
        }
    }
    inside
}
fn in_face(p: [f64; 2], f: &FlatFace) -> bool {
    f.loops.iter().filter(|lp| in_loop(p, lp)).count() % 2 == 1
}

struct Grid {
    coordinates: [Vec<f64>; 3],
    centers: [Vec<f64>; 3],
    dimensions: [usize; 3],
}
impl Grid {
    fn new(a: &Body, b: &Body, tolerance: f64) -> Result<Self> {
        let mut coordinates: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        for p in a.vertices.iter().chain(&b.vertices) {
            for i in 0..3 {
                coordinates[i].push(xyz(*p)[i]);
            }
        }
        let mut centers: [Vec<f64>; 3] = std::array::from_fn(|_| Vec::new());
        for i in 0..3 {
            coordinates[i].sort_by(f64::total_cmp);
            coordinates[i].dedup_by(|a, b| *a == *b);
            if coordinates[i].len() > MAX_COORDINATES {
                return Err(unsupported("coordinate-plane limit exceeded"));
            }
            if coordinates[i].len() < 2 {
                return Err(invalid("body has no enclosed extent"));
            }
            for pair in coordinates[i].windows(2) {
                let d = pair[1] - pair[0];
                if !d.is_finite() || d < tolerance {
                    return Err(unsupported(
                        "distinct planes below supported tolerance separation",
                    ));
                }
                let mid = pair[0] + d * 0.5;
                if mid <= pair[0] || mid >= pair[1] {
                    return Err(unsupported("no representable interval interior"));
                }
                centers[i].push(mid);
            }
        }
        let dimensions = std::array::from_fn(|i| centers[i].len());
        if dimensions.iter().product::<usize>() > MAX_CELLS {
            return Err(unsupported("rectilinear cell limit exceeded"));
        }
        Ok(Self {
            coordinates,
            centers,
            dimensions,
        })
    }
    fn index(&self, p: [usize; 3]) -> usize {
        (p[0] * self.dimensions[1] + p[1]) * self.dimensions[2] + p[2]
    }
    fn membership(&self, faces: &[FlatFace]) -> Result<Vec<bool>> {
        let mut cells = vec![false; self.dimensions.iter().product()];
        let mut tests = 0;
        let mut groups: [Vec<Vec<&FlatFace>>; 3] =
            std::array::from_fn(|i| vec![Vec::new(); self.coordinates[i].len()]);
        for f in faces {
            let at = self.coordinates[f.axis]
                .iter()
                .position(|v| *v == f.coordinate)
                .ok_or_else(|| invalid("plane absent from vertex coordinates"))?;
            groups[f.axis][at].push(f);
        }
        for y in 0..self.dimensions[1] {
            for z in 0..self.dimensions[2] {
                let mut inside = false;
                for x in 0..self.coordinates[0].len() {
                    let crossing = coverage(
                        [self.centers[1][y], self.centers[2][z]],
                        &groups[0][x],
                        &mut tests,
                    )?;
                    match crossing {
                        -1 if !inside => inside = true,
                        1 if inside => inside = false,
                        0 => {}
                        _ => return Err(invalid("overlapping or inconsistently oriented shell")),
                    }
                    if x < self.dimensions[0] {
                        cells[self.index([x, y, z])] = inside;
                    }
                }
                if inside {
                    return Err(invalid("body scan did not close"));
                }
            }
        }
        // Independent surface coverage in all directions verifies that the
        // inferred material exactly matches every declared operand boundary.
        for axis in 0..3 {
            let u = (axis + 1) % 3;
            let v = (axis + 2) % 3;
            for at in 0..self.coordinates[axis].len() {
                for i in 0..self.dimensions[u] {
                    for j in 0..self.dimensions[v] {
                        let mut p = [0; 3];
                        p[u] = i;
                        p[v] = j;
                        let before = if at == 0 {
                            false
                        } else {
                            p[axis] = at - 1;
                            cells[self.index(p)]
                        };
                        let after = if at == self.dimensions[axis] {
                            false
                        } else {
                            p[axis] = at;
                            cells[self.index(p)]
                        };
                        let actual = coverage(
                            [self.centers[u][i], self.centers[v][j]],
                            &groups[axis][at],
                            &mut tests,
                        )?;
                        if actual != i8::from(before) - i8::from(after) {
                            return Err(invalid("face coverage disagrees with enclosed material"));
                        }
                    }
                }
            }
        }
        Ok(cells)
    }
    fn require_connected_material(&self, cells: &[bool]) -> Result<()> {
        let Some(first) = cells.iter().position(|v| *v) else {
            return Ok(());
        };
        let z = first % self.dimensions[2];
        let t = first / self.dimensions[2];
        let mut stack = vec![[t / self.dimensions[1], t % self.dimensions[1], z]];
        let mut seen = vec![false; cells.len()];
        seen[first] = true;
        while let Some(p) = stack.pop() {
            for axis in 0..3 {
                for positive in [false, true] {
                    let mut n = p;
                    if positive {
                        if n[axis] + 1 >= self.dimensions[axis] {
                            continue;
                        }
                        n[axis] += 1;
                    } else {
                        if n[axis] == 0 {
                            continue;
                        }
                        n[axis] -= 1;
                    }
                    let i = self.index(n);
                    if cells[i] && !seen[i] {
                        seen[i] = true;
                        stack.push(n);
                    }
                }
            }
        }
        if cells.iter().zip(seen).any(|(cell, seen)| *cell && !seen) {
            return Err(unsupported(
                "multiple disconnected result solids are not supported yet",
            ));
        }
        Ok(())
    }
    fn make_body(&self, cells: &[bool]) -> Body {
        let mut body = Body::default();
        let mut vertices = BTreeMap::<[usize; 3], usize>::new();
        let mut edges = BTreeMap::<(usize, usize), usize>::new();
        for x in 0..self.dimensions[0] {
            for y in 0..self.dimensions[1] {
                for z in 0..self.dimensions[2] {
                    let cell = [x, y, z];
                    if !cells[self.index(cell)] {
                        continue;
                    }
                    for axis in 0..3 {
                        for positive in [false, true] {
                            let mut next = cell;
                            let exists = if positive {
                                next[axis] += 1;
                                next[axis] < self.dimensions[axis]
                            } else if next[axis] > 0 {
                                next[axis] -= 1;
                                true
                            } else {
                                false
                            };
                            if exists && cells[self.index(next)] {
                                continue;
                            }
                            let u = (axis + 1) % 3;
                            let v = (axis + 2) % 3;
                            let mut corner = cell;
                            if positive {
                                corner[axis] += 1;
                            }
                            let mut corners = [corner; 4];
                            corners[1][u] += 1;
                            corners[2][u] += 1;
                            corners[2][v] += 1;
                            corners[3][v] += 1;
                            if !positive {
                                corners.reverse();
                            }
                            let indices: Vec<usize> = corners
                                .iter()
                                .map(|p| {
                                    *vertices.entry(*p).or_insert_with(|| {
                                        let i = body.vertices.len();
                                        body.vertices.push(vec(std::array::from_fn(|a| {
                                            self.coordinates[a][p[a]]
                                        })));
                                        i
                                    })
                                })
                                .collect();
                            let mut lp = Vec::new();
                            for i in 0..4 {
                                let a = indices[i];
                                let b = indices[(i + 1) % 4];
                                let key = (a.min(b), a.max(b));
                                let edge = *edges.entry(key).or_insert_with(|| {
                                    let i = body.edges.len();
                                    body.edges.push(Edge {
                                        curve: Curve::Line {
                                            a: body.vertices[key.0],
                                            b: body.vertices[key.1],
                                        },
                                        start: key.0,
                                        end: key.1,
                                    });
                                    i
                                });
                                lp.push((edge, a == key.0));
                            }
                            let mut n = [0.; 3];
                            n[axis] = if positive { 1. } else { -1. };
                            let mut reference = [0.; 3];
                            reference[u] = 1.;
                            body.faces.push(Face {
                                surface: Surface::Plane {
                                    p: body.vertices[indices[0]],
                                    n: vec(n),
                                    x: vec(reference),
                                },
                                sense: true,
                                loops: vec![lp],
                            });
                        }
                    }
                }
            }
        }
        body
    }
}
fn coverage(p: [f64; 2], faces: &[&FlatFace], tests: &mut usize) -> Result<i8> {
    let mut result = 0;
    for face in faces {
        *tests += 1;
        if *tests > MAX_CHECKS {
            return Err(unsupported("boundary check budget exceeded"));
        }
        if in_face(p, face) {
            if result != 0 {
                return Err(invalid("overlapping coplanar face interiors"));
            }
            result = face.sign;
        }
    }
    Ok(result)
}

#[cfg(test)]
fn cut_rectilinear(a: &Body, b: &Body) -> Result<Body> {
    cut(a, b, 1e-7)
}
#[cfg(test)]
#[path = "tests.rs"]
mod tests;
