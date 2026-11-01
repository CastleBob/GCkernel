//! Offset and shell solids (Claude 2026-11-01; Bob: "continue with ... shell and offset solids"). Exact construction from
//! the public mathematics of parallel bodies (README.md): every face of the result is an exact offset of a source face
//! (a parallel plane, a coaxial cylinder or cone, a concentric sphere, a torus with its tube radius changed) or, in a
//! rounded outward offset, the exact rolling-ball surface about a source edge or vertex (cylinder, torus, sphere).
//!
//! Domains, each proved from the body's topology and geometry and refused by name otherwise:
//! - straight prisms of lines and arcs (plates, plates with holes, tube, angle, channel and beam profiles): prism.rs;
//! - full solids of revolution of lines and arcs (cylinders, cones, spheres, tori, revolved sections): revolved.rs.
//!
//! OFFSET moves every face by a signed distance (positive outward): `Rounded` is the true offset (the set of points
//! within the distance of the solid, or at least that far inside it), `Sharp` extends the moved faces to meet instead
//! of rounding. SHELL hollows the solid to a wall thickness (positive inward, negative outward): with no face removed
//! the cavity is a void, an inner shell of the same body; removing faces opens the cavity through them.
//!
//! Nothing partial is returned: a face that would vanish, a cavity that would split or merge, an unsupported surface or
//! domain is refused by name. Every result passes `topology::check` and its exact volume (closed form from the
//! construction) agrees with `mass::measure` of the built body, or the operation reports an internal fault.

pub mod planar;
pub mod prism;
pub mod revolved;

use crate::neutral::{write_body, Body, Reader};
use crate::{mass, topology, KernelError};

/// What a result face was made from: a source face (offset or kept), or the source edge or vertex a rounding face is
/// the rolling-ball surface about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    Face(usize),
    Edge(usize),
    Vertex(usize),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Rounded,
    Sharp,
}

/// A result: the body, its exact volume, and per face its source.
#[derive(Clone, Debug)]
pub struct Offset {
    pub body: Body,
    pub volume: f64,
    pub sources: Vec<Source>,
}

fn unsupported<T>(m: String) -> Result<T, KernelError> { Err(KernelError::Unsupported(m)) }

/// The body turned inside out: every face's sense flipped and its loops run backwards (a void's boundary).
pub fn flipped(body: &Body) -> Body {
    let mut out = body.clone();
    for f in &mut out.faces {
        f.sense = !f.sense;
        for lp in &mut f.loops {
            lp.reverse();
            for s in lp.iter_mut() { s.1 = !s.1; }
        }
    }
    out
}

/// `extra` added to `body` (its vertices, edges and faces renumbered after the existing ones).
pub fn append(body: &mut Body, extra: &Body) {
    let (vo, eo) = (body.vertices.len(), body.edges.len());
    body.vertices.extend_from_slice(&extra.vertices);
    body.edges.extend(extra.edges.iter().map(|e| crate::neutral::Edge { curve: e.curve.clone(), start: e.start + vo, end: e.end + vo }));
    body.faces.extend(extra.faces.iter().map(|f| crate::neutral::Face {
        surface: f.surface.clone(), sense: f.sense, loops: f.loops.iter().map(|lp| lp.iter().map(|&(e, d)| (e + eo, d)).collect()).collect(),
    }));
}

/// The domain a body was proved to be in.
enum Domain {
    Prism(prism::Prism),
    Revolved(revolved::Revolved),
}

fn domain(body: &Body, caps: &[usize], tol: f64) -> Result<Domain, KernelError> {
    let why_prism = match prism::recognize(body, caps, tol) {
        Ok(p) => return Ok(Domain::Prism(p)),
        Err(KernelError::Unsupported(m)) => m,
        Err(e) => return Err(e),
    };
    let why_revolved = match revolved::recognize(body, tol) {
        Ok(r) => return Ok(Domain::Revolved(r)),
        Err(KernelError::Unsupported(m)) => m,
        Err(e) => return Err(e),
    };
    unsupported(format!(
        "offsets and shells are built for straight prisms and full solids of revolution of lines and arcs; this body is {why_prism}, \
         and {why_revolved}"))
}

/// Built results are checked as every kernel operation's are: a closed, consistent solid whose measured volume is the
/// construction's closed form.
fn finish(body: Body, volume: f64, sources: Vec<Source>, tol: f64, what: &str) -> Result<Offset, KernelError> {
    topology::check(&body, tol).map_err(|e| KernelError::Internal(format!("{what} built an invalid solid: {e}")))?;
    let measured = mass::measure(&body, tol).map_err(|e| KernelError::Internal(format!("{what} built a solid that cannot be measured: {e}")))?;
    if (measured.volume - volume).abs() > 1e-9 * volume.abs().max(1e-300) {
        return Err(KernelError::Internal(format!("{what}: the built solid measures {} but its construction gives {volume}", measured.volume)));
    }
    if sources.len() != body.faces.len() {
        return Err(KernelError::Internal(format!("{what}: {} face sources for {} faces", sources.len(), body.faces.len())));
    }
    Ok(Offset { body, volume, sources })
}

/// Every face moved by `distance` (positive outward).
pub fn offset(body: &Body, distance: f64, mode: Mode, tol: f64) -> Result<Offset, KernelError> {
    if !(distance.is_finite() && distance != 0.0) { return Err(KernelError::InvalidInput("an offset needs a finite, non-zero distance".into())); }
    if distance.abs() <= tol { return Err(KernelError::InvalidInput(format!("an offset of {distance} is within the tolerance {tol}"))); }
    topology::check(body, tol)?;
    let sharp = mode == Mode::Sharp;
    match domain(body, &[], tol)? {
        Domain::Prism(p) => {
            let (b, roles, volume) = prism::offset(&p, distance, sharp, tol)?;
            let sources = roles.iter().map(|r| p.source(r)).collect();
            finish(b, volume, sources, tol, "offset")
        }
        Domain::Revolved(r) => {
            let (b, sources, volume) = revolved::offset(body, &r, distance, sharp, tol)?;
            finish(b, volume, sources, tol, "offset")
        }
    }
}

/// The solid hollowed to wall `thickness` (positive inward, negative outward), open through the `removed` faces.
pub fn shell(body: &Body, thickness: f64, removed: &[usize], tol: f64) -> Result<Offset, KernelError> {
    if !(thickness.is_finite() && thickness != 0.0) { return Err(KernelError::InvalidInput("a shell needs a finite, non-zero thickness".into())); }
    if thickness.abs() <= tol { return Err(KernelError::InvalidInput(format!("a shell wall of {thickness} is within the tolerance {tol}"))); }
    topology::check(body, tol)?;
    for (i, &f) in removed.iter().enumerate() {
        if f >= body.faces.len() { return Err(KernelError::InvalidInput(format!("removed face {f} is not one of the body's {} faces", body.faces.len()))); }
        if removed[..i].contains(&f) { return Err(KernelError::InvalidInput(format!("face {f} is removed twice"))); }
    }
    if removed.is_empty() {
        // A closed shell: the body and its offset, the inner one a void.
        let measured = mass::measure(body, tol)?;
        if !measured.exact { return unsupported("a closed shell needs a body measured in closed form".into()); }
        let moved = offset(body, -thickness, Mode::Rounded, tol)?;
        let own: Vec<Source> = (0..body.faces.len()).map(Source::Face).collect();
        let (b, sources, volume) = if thickness > 0.0 {
            let mut b = body.clone();
            append(&mut b, &flipped(&moved.body));
            let mut s = own;
            s.extend(moved.sources);
            (b, s, measured.volume - moved.volume)
        } else {
            let mut b = moved.body.clone();
            append(&mut b, &flipped(body));
            let mut s = moved.sources;
            s.extend(own);
            (b, s, moved.volume - measured.volume)
        };
        return finish(b, volume, sources, tol, "shell");
    }
    match domain(body, removed, tol)? {
        Domain::Prism(p) => {
            let (b, roles, volume) = prism::shell_open(&p, thickness, removed.contains(&p.top), removed.contains(&p.bottom), tol)?;
            let sources = roles.iter().map(|r| p.source(r)).collect();
            finish(b, volume, sources, tol, "shell")
        }
        Domain::Revolved(r) => {
            let (b, sources, volume) = revolved::shell_open(body, &r, thickness, removed, tol)?;
            finish(b, volume, sources, tol, "shell")
        }
    }
}

/// The reply shared by ops 20 and 21: the body, its exact volume, the face count and per face [kind, index]
/// (kind 0 a source face, 1 a source edge, 2 a source vertex).
fn reply(result: &Offset) -> Vec<f64> {
    let mut out = write_body(&result.body);
    out.push(result.volume);
    out.push(result.sources.len() as f64);
    for s in &result.sources {
        let (kind, index) = match *s { Source::Face(i) => (0.0, i), Source::Edge(i) => (1.0, i), Source::Vertex(i) => (2.0, i) };
        out.extend_from_slice(&[kind, index as f64]);
    }
    out
}

/// Op 20 SHELL: [thickness, k, k face indices, body] -> `reply`.
pub fn shell_op(data: &[f64], tol: f64) -> Result<Vec<f64>, KernelError> {
    let mut reader = Reader::with_tolerance(data, tol);
    let thickness = reader.next()?;
    let k = reader.count()?;
    let mut removed = Vec::with_capacity(k);
    for _ in 0..k {
        let v = reader.next()?;
        if v < 0.0 || v.fract() != 0.0 { return Err(KernelError::InvalidInput(format!("removed face index {v} is not a whole number"))); }
        removed.push(v as usize);
    }
    let body = reader.body()?;
    reader.finish()?;
    Ok(reply(&shell(&body, thickness, &removed, tol)?))
}

/// Op 21 OFFSET: [distance, mode (0 rounded, 1 sharp), body] -> `reply`.
pub fn offset_op(data: &[f64], tol: f64) -> Result<Vec<f64>, KernelError> {
    let mut reader = Reader::with_tolerance(data, tol);
    let distance = reader.next()?;
    let mode = match reader.next()? {
        m if m == 0.0 => Mode::Rounded,
        m if m == 1.0 => Mode::Sharp,
        m => return Err(KernelError::InvalidInput(format!("offset mode must be 0 (rounded) or 1 (sharp), not {m}"))),
    };
    let body = reader.body()?;
    reader.finish()?;
    Ok(reply(&offset(&body, distance, mode, tol)?))
}
