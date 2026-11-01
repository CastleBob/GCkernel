//! Radial bores through the walls of revolved stock with conical faces (Claude 2026-10-31, tapered Pole stock, solid
//! and hollow).
//!
//! Domain, decided from the geometry before anything is built: a full round cutter (one cylinder, two caps) whose axis
//! is square to a solid of revolution's axis; the walls it crosses are each ONE face, a cone or cylinder about that
//! axis, whose axial extent holds the hole's whole axial window [z - rho, z + rho] strictly inside its end rings. The
//! cutter starts outside the outer wall and ends in the bore (one wall), past the far side (through), or, in solid
//! stock, inside the material past every entry point (blind). Each crossing is one closed loop on its wall: on a cone
//! the exact cylinder/cone intersection, `Curve::ConePair` (wire type 8: per angle round the hole the quadratic for
//! the height along its axis, discriminant proved positive all the way round), on a cylinder a CylinderPair. The
//! cutter's cylinder joins the loops (a blind hole: its entry loop to its bottom disk). Nothing else changes.
//!
//! Clearances (independent, elementary, d square to the axis b): every point of a loop lies at axial height within
//! rho of the hole's axis (h d has no part along b), and at an angle about the axis within asin((|w| + rho) / R_min)
//! of the hole's direction (w its axis's offset across both, R_min the wall's least radius over the window), so
//! rings, seams and earlier holes on a wall are cleared by those windows; every other face's edges must be clear of
//! the whole cutter's box (rings and lines by their extrema, pair edges by the bound below). Anything else is refused
//! by name; nothing is retried.
use crate::geom::{cone_pair_disc_coefficients, pair_disc_coefficients, pair_disc_lower_bound, Curve, Surface, Vec3};
use crate::neutral::{Body, Edge, Face};
use crate::KernelError;
use std::f64::consts::{PI, TAU};

type Result<T> = std::result::Result<T, KernelError>;

fn unsupported(message: &str) -> KernelError {
    KernelError::Unsupported(format!("bore through a conical wall: {message}"))
}

fn perp(v: Vec3, b: Vec3) -> Vec3 { v.sub(b.scale(v.dot(b))) }

/// A wall: its face, and its radius as r0 + k (z - z0) along the axis (k = 0 on a cylinder).
#[derive(Clone, Copy)]
struct Wall {
    face: usize,
    cone: bool,
    /// The surface's own frame: point, axis, radius there, slope.
    sp: Vec3,
    sa: Vec3,
    sr: f64,
    sk: f64,
}

impl Wall {
    /// The radius at axial height z (from the body axis point p along b).
    fn radius(&self, p: Vec3, b: Vec3, z: f64) -> f64 {
        let s = p.add(b.scale(z)).sub(self.sp).dot(self.sa);
        self.sr + self.sk * s
    }
    fn range(&self, p: Vec3, b: Vec3, lo: f64, hi: f64) -> (f64, f64) {
        let (a, c) = (self.radius(p, b, lo), self.radius(p, b, hi));
        (a.min(c), a.max(c))
    }
    /// The loop where the cutter (start s, direction d, radius rho; seam frame x) crosses this wall.
    fn loop_curve(&self, s: Vec3, d: Vec3, x: Vec3, rho: f64, branch: f64) -> Curve {
        if self.cone {
            Curve::ConePair { c: s, a: d, x, rho, p: self.sp, b: self.sa, r: self.sr, k: self.sk, branch, t0: 0., t1: TAU }
        } else {
            Curve::CylinderPair { c: s, a: d, x, rho, p: self.sp, b: self.sa, r: self.sr, branch, t0: 0., t1: TAU }
        }
    }
    fn normal(&self, q: Vec3, sense: bool) -> Result<Vec3> {
        let w = perp(q.sub(self.sp), self.sa).unit()?;
        let n = w.sub(self.sa.scale(self.sk)).unit()?;
        Ok(n.scale(if sense { 1. } else { -1. }))
    }
}

/// The wall surfaces about the axis line (p, b): (face, wall) for every cone or cylinder face on it, with its sense.
fn walls(body: &Body, p: Vec3, b: Vec3, tol: f64) -> Vec<(Wall, bool)> {
    body.faces.iter().enumerate().filter_map(|(i, f)| {
        let (sp, sa, sr, sk, cone) = match f.surface {
            Surface::Cone { p, a, r, k, .. } => (p, a, r, k, true),
            Surface::Cylinder { p, a, r, .. } => (p, a, r, 0., false),
            _ => return None,
        };
        (sa.cross(b).length() <= 1e-9 && perp(sp.sub(p), b).length() <= tol)
            .then_some((Wall { face: i, cone, sp, sa, sr, sk }, f.sense))
    }).collect()
}

/// The axial extent of a wall face (its edges' heights along b from p): rings give their level, lines their ends.
fn extent(body: &Body, face: usize, p: Vec3, b: Vec3) -> Option<(f64, f64)> {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for (e, _) in body.faces[face].loops.iter().flatten() {
        let z = |q: Vec3| q.sub(p).dot(b);
        match body.edges[*e].curve {
            Curve::Arc { c, n, .. } if n.cross(b).length() <= 1e-9 => { lo = lo.min(z(c)); hi = hi.max(z(c)); }
            Curve::Line { a, b: e2 } => { lo = lo.min(z(a).min(z(e2))); hi = hi.max(z(a).max(z(e2))); }
            Curve::CylinderPair { .. } | Curve::ConePair { .. } => {}
            _ => return None,
        }
    }
    lo.is_finite().then_some((lo, hi))
}

/// The cutter and the walls it crosses, after every admission and clearance check.
struct Plan {
    s: Vec3,
    d: Vec3,
    rho: f64,
    len: f64,
    b: Vec3,
    /// (wall, wall face sense, branch) outside in.
    crossings: Vec<(Wall, bool, f64)>,
    /// A blind hole's bottom disk centre.
    blind: Option<Vec3>,
}

fn plan(body: &Body, tool: &Body, tol: f64) -> Result<Plan> {
    let (c0, c1, rho) = super::radial::cutter(tool, tol).ok_or_else(|| unsupported("the tool is not a full round cutter"))?;
    let axis = c1.sub(c0).unit()?;
    // The revolved body's axis: its cone faces' (all one axis line).
    let cones: Vec<(Vec3, Vec3)> = body.faces.iter().filter_map(|f| match f.surface {
        Surface::Cone { p, a, .. } => Some((p, a)),
        _ => None,
    }).collect();
    let Some(&(p, b)) = cones.first() else { return Err(unsupported("the body has no conical wall")) };
    if cones.iter().any(|(q, a)| a.cross(b).length() > 1e-9 || perp(q.sub(p), b).length() > tol) {
        return Err(unsupported("the body's conical faces are not on one axis"));
    }
    if axis.dot(b).abs() > 1e-9 {
        return Err(unsupported("the hole is not square to the stock's axis"));
    }
    // Outside in.
    let dist = |q: Vec3| perp(q.sub(p), b).length();
    let (s, e) = if dist(c0) >= dist(c1) { (c0, c1) } else { (c1, c0) };
    let d = e.sub(s).unit()?;
    let len = e.sub(s).length();
    let z = s.sub(p).dot(b);
    let (zlo, zhi) = (z - rho, z + rho);
    let offset = s.sub(p).dot(d.cross(b)).abs();
    let (sd, ed) = (s.sub(p).dot(d), e.sub(p).dot(d));
    let candidates = walls(body, p, b, tol);
    let holding = |w: &Wall| extent(body, w.face, p, b).is_some_and(|(lo, hi)| lo < zlo - tol && hi > zhi + tol);
    let outer: Vec<&(Wall, bool)> = candidates.iter().filter(|(w, sense)| *sense && holding(w)).collect();
    let inner: Vec<&(Wall, bool)> = candidates.iter().filter(|(w, sense)| !*sense && holding(w)).collect();
    let [&(wo, _)] = outer[..] else { return Err(unsupported("the hole's axial window is not inside one outer wall face")) };
    if inner.len() > 1 {
        return Err(unsupported("the hole's axial window meets more than one inner wall"));
    }
    let (ro_min, ro_max) = wo.range(p, b, zlo, zhi);
    if !(-sd > ro_max + tol && sd < 0.) {
        return Err(unsupported("the cutter does not start outside the outer wall"));
    }
    let mut crossings = vec![(wo, true, -1.)];
    let mut blind = None;
    let reach = |r: f64| (ed * ed + (offset + rho).powi(2)).sqrt() < r - tol;
    match inner.first() {
        Some(&&(wi, _)) => {
            let (ri_min, _) = wi.range(p, b, zlo, zhi);
            if offset + rho >= ri_min - tol {
                return Err(unsupported("the hole is not clear of tangency with the bore"));
            }
            crossings.push((wi, false, -1.));
            if ed > ro_max + tol {
                crossings.extend([(wi, false, 1.), (wo, true, 1.)]);
            } else if !reach(ri_min) {
                return Err(unsupported("the cutter ends inside the wall (neither in the bore nor past the far side)"));
            }
        }
        None => {
            if offset + rho >= ro_min - tol {
                return Err(unsupported("the hole is not clear of tangency with the outer wall"));
            }
            if ed > ro_max + tol {
                crossings.push((wo, true, 1.));
            } else {
                // Blind in solid stock: the bottom disk inside the material, deeper than every entry point (each at
                // -sqrt(R^2 - w^2) <= -sqrt(R_min^2 - (offset + rho)^2) along d from the axis).
                let deepest_entry = -(ro_min * ro_min - (offset + rho).powi(2)).sqrt();
                if !(reach(ro_min) && ed > deepest_entry + tol) {
                    return Err(unsupported("the blind hole's bottom is not inside the solid past its entry"));
                }
                blind = Some(e);
            }
        }
    }
    // Discriminants all the way round.
    let x = b;
    for &(w, _, _) in &crossings {
        let q = if w.cone {
            cone_pair_disc_coefficients(s, d, x, rho, w.sp, w.sa, w.sr, w.sk)
        } else {
            pair_disc_coefficients(s, d, x, rho, w.sp, w.sa, w.sr)
        };
        if !(pair_disc_lower_bound(&q, 0., TAU, 4096) > tol * tol) {
            return Err(unsupported("a loop is not proved clear of tangency"));
        }
    }
    // Each crossed wall: its rings, seams and earlier holes clear of the new loops (axial window or angle).
    let e1 = b.perpendicular();
    let e2 = b.cross(e1);
    let angle_of = |v: Vec3| v.dot(e2).atan2(v.dot(e1));
    let interval = |centre: Vec3, half: f64| { let a = angle_of(centre); (a - half, a + half) };
    let meet = |i: (f64, f64), j: (f64, f64), margin: f64| {
        let dm = ((i.0 + i.1) / 2. - (j.0 + j.1) / 2. + PI).rem_euclid(TAU) - PI;
        dm.abs() <= (i.1 - i.0) / 2. + (j.1 - j.0) / 2. + margin
    };
    for &(w, _, _) in &crossings {
        let (rmin, _) = w.range(p, b, zlo, zhi);
        let reach = (offset + rho) / rmin;
        let new: Vec<(f64, f64)> = crossings.iter().filter(|c| c.0.face == w.face)
            .map(|c| interval(d.scale(c.2), reach.asin())).collect();
        for (e, _) in body.faces[w.face].loops.iter().flatten() {
            let zof = |q: Vec3| q.sub(p).dot(b);
            match body.edges[*e].curve {
                Curve::Arc { c, .. } => {
                    if (zof(c) - z).abs() <= rho + tol {
                        return Err(unsupported("the hole meets an end of its wall"));
                    }
                }
                Curve::Line { a, b: q } => {
                    let (za, zq) = (zof(a), zof(q));
                    let axial = za.min(zq) > zhi + tol || za.max(zq) < zlo - tol;
                    let theta = angle_of(perp(a.sub(p), b));
                    if !axial && new.iter().any(|i| meet(*i, (theta, theta), tol / rmin)) {
                        return Err(unsupported("the hole meets its wall's seam"));
                    }
                }
                Curve::CylinderPair { c, a, rho: hr, branch, .. } | Curve::ConePair { c, a, rho: hr, branch, .. } => {
                    if a.dot(b).abs() > 1e-9 {
                        return Err(unsupported("an earlier hole on the wall is not square to the axis"));
                    }
                    let hz = c.sub(p).dot(b);
                    let axial = hz - hr > zhi + tol || hz + hr < zlo - tol;
                    let hoff = c.sub(p).dot(a.cross(b)).abs();
                    let (hmin, _) = w.range(p, b, hz - hr, hz + hr);
                    if !((hoff + hr) / hmin < 1.) {
                        return Err(unsupported("an earlier hole on the wall has no angular bound"));
                    }
                    let old = interval(a.scale(branch), ((hoff + hr) / hmin).asin());
                    if !axial && new.iter().any(|i| meet(*i, old, tol / rmin)) {
                        return Err(unsupported("the hole meets an earlier hole on its wall"));
                    }
                }
                _ => return Err(unsupported("a wall edge is neither a ring, a seam nor a hole")),
            }
        }
    }
    // Every other face clear of the whole cutter.
    let cutter_box = {
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for e in &tool.edges {
            let (l, h) = super::passage::curve_bounds(&e.curve);
            for k in 0..3 { lo[k] = lo[k].min(l[k]); hi[k] = hi[k].max(h[k]); }
        }
        (lo, hi)
    };
    let crossed: Vec<usize> = crossings.iter().map(|c| c.0.face).collect();
    for (i, f) in body.faces.iter().enumerate() {
        if crossed.contains(&i) {
            continue;
        }
        for (e, _) in f.loops.iter().flatten() {
            let bx = edge_box(&body.edges[*e].curve, &candidates).ok_or_else(|| unsupported("an edge has no bounding box"))?;
            if !(0..3).any(|k| bx.1[k] + tol < cutter_box.0[k] || cutter_box.1[k] + tol < bx.0[k]) {
                return Err(unsupported("the hole meets another face of the body"));
            }
        }
    }
    Ok(Plan { s, d, rho, len, b, crossings, blind })
}

/// A rigorous box: lines and arcs by their extrema; a pair edge square to its wall's axis: across both axes within
/// rho of its hole's axis, along the axis b within rho, and along the hole on its wall at sqrt(R^2 - w^2) for the
/// wall's radius range and |w| at most the offset plus rho.
fn edge_box(curve: &Curve, walls: &[(Wall, bool)]) -> Option<([f64; 3], [f64; 3])> {
    match *curve {
        Curve::Line { .. } | Curve::Arc { .. } => Some(super::passage::curve_bounds(curve)),
        Curve::CylinderPair { c, a, rho, p, b, r, branch, .. } => cross_box(c, a, rho, p, b, branch, (r, r)),
        Curve::ConePair { c, a, rho, p, b, r, k, branch, .. } => {
            let _ = walls;
            let z = c.sub(p).dot(b);
            let (r0, r1) = (r + k * (z - rho), r + k * (z + rho));
            cross_box(c, a, rho, p, b, branch, (r0.min(r1), r0.max(r1)))
        }
        _ => None,
    }
}

#[allow(clippy::too_many_arguments)]
fn cross_box(c: Vec3, a: Vec3, rho: f64, p: Vec3, b: Vec3, branch: f64, (rmin, rmax): (f64, f64)) -> Option<([f64; 3], [f64; 3])> {
    if a.dot(b).abs() > 1e-9 {
        return None;
    }
    let w = b.cross(a);
    let (cb, cw) = (c.sub(p).dot(b), c.sub(p).dot(w));
    let (w0, w1) = (cw - rho, cw + rho);
    let m = w0.abs().max(w1.abs());
    if !(m < rmin) {
        return None;
    }
    let near = (rmin * rmin - m * m).sqrt();
    let along = if branch > 0. { [near, rmax] } else { [-rmax, -near] };
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for sb in [cb - rho, cb + rho] {
        for sw in [w0, w1] {
            for sa in along {
                let q = p.add(b.scale(sb)).add(w.scale(sw)).add(a.scale(sa));
                for (k, v) in [q.x, q.y, q.z].into_iter().enumerate() {
                    lo[k] = lo[k].min(v);
                    hi[k] = hi[k].max(v);
                }
            }
        }
    }
    Some((lo, hi))
}

/// Admission (when the body has conical faces): every condition holds, or the named reason it does not.
pub(super) fn admit(body: &Body, tool: &Body, tol: f64) -> Result<()> { plan(body, tool, tol).map(|_| ()) }

pub(super) fn cut(body: &Body, tool: &Body, tol: f64) -> Result<Body> {
    crate::topology::check(body, tol)?;
    crate::topology::check(tool, tol)?;
    super::radial::planar_faces_outward(body)?;
    super::radial::planar_faces_outward(tool)?;
    let plan = plan(body, tool, tol)?;
    let (s, d, rho, b) = (plan.s, plan.d, plan.rho, plan.b);
    let mut out = body.clone();
    let mut flags = vec![];
    for &(w, sense, branch) in &plan.crossings {
        let curve = w.loop_curve(s, d, b, rho, branch);
        let (q, tangent) = curve.eval(0.);
        let v = out.vertices.len();
        out.vertices.push(q);
        let e = out.edges.len();
        out.edges.push(Edge { curve, start: v, end: v });
        // A hole runs clockwise about the wall's outward normal: the hole (toward the cutter's axis) on the right.
        let normal = w.normal(q, sense)?;
        let foot = s.add(d.scale(q.sub(s).dot(d)));
        let inward = perp(foot.sub(q), normal);
        let forward = tangent.cross(normal).dot(inward) > 0.;
        out.faces[w.face].loops.push(vec![(e, forward)]);
        flags.push((e, forward, v));
    }
    let wall = |out: &mut Body, (e_in, f_in, v_in): (usize, bool, usize), (e_out, f_out, v_out): (usize, bool, usize)| {
        let seam = out.edges.len();
        out.edges.push(Edge { curve: Curve::Line { a: out.vertices[v_in], b: out.vertices[v_out] }, start: v_in, end: v_out });
        out.faces.push(Face {
            surface: Surface::Cylinder { p: s, a: d, x: b, r: rho },
            sense: false,
            loops: vec![vec![(e_in, !f_in), (seam, true), (e_out, f_out), (seam, false)]],
        });
    };
    match plan.blind {
        None => {
            for pair in flags.chunks_exact(2) {
                wall(&mut out, pair[0], (pair[1].0, !pair[1].1, pair[1].2));
            }
        }
        Some(e) => {
            // The bottom disk at the cutter's end (facing back up the hole) and the hole's wall down to it.
            let start = e.add(b.scale(rho));
            let v = out.vertices.len();
            out.vertices.push(start);
            let ring = out.edges.len();
            out.edges.push(Edge { curve: Curve::Arc { c: e, n: d, x: b, r: rho, t0: 0., t1: TAU }, start: v, end: v });
            // The entry loop runs counter-clockwise about d in the cutter's face exactly when !forward; the bottom ring
            // the other way round there.
            let (e_in, f_in, v_in) = flags[0];
            wall(&mut out, (e_in, f_in, v_in), (ring, f_in, v));
            out.faces.push(Face { surface: Surface::Plane { p: e, n: d.scale(-1.), x: b }, sense: true, loops: vec![vec![(ring, !f_in)]] });
        }
    }
    let _ = plan.len;
    crate::topology::check(&out, tol)?;
    Ok(out)
}
