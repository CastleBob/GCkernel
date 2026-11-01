# Fillets and chamfers (blend)

Owner: Claude (2026-11-01, Bob: "continue with ... the fillets and the chamfers"). Clean room: written from the
public mathematics below; no third-party kernel source was read or translated. No new dependencies.

Entry: `blend(body, items, tolerance) -> Result<Blended, KernelError>`; worker op 19 (`protocol.rs`). Result faces
0..F-1 are the input faces in their order; then one blend face per edge (input edge order); then corner patches.

## Supported exactly

| Edge | Chamfer | Fillet |
|---|---|---|
| straight, between two flat faces (convex or concave, any opening angle) | plane | circular cylinder |
| circle where a flat face meets a cylinder at right angles (hole rim, bar/pipe end, shoulder) | cone | torus |

Chamfer forms: equal distance; two distances (the first on a named face); distance and angle (from a named face).

Ends, at vertices where three faces meet:
- one blend there: it runs out on the flat end face (section: line, circle, or ellipse for a slanted face);
- two: a mitre (chamfer planes meet in a line; equal fillets in an ellipse) or, along a rim, the torus meridian or
  cone generator through the seam;
- three: chamfers get the triangular plane through the face corners; equal fillets get a sphere patch when one face
  is square to the other two (box corners, any extruded plate's top corners).

## Refused by name

Vertices where other than three faces meet (including the coplanar seams the rectilinear cut leaves); blends ending on
curved faces; fillet + chamfer, convex + concave, or unequal fillets meeting at a vertex; mitres whose sides reach the
third edge at different points; three-fillet corners without a square face; rims blended part way round or unalike;
planes not square to the cylinder; edges on other surfaces; any blend that consumes a face or edge, crosses another
loop, reaches the axis, runs past a cylinder face's other boundary, or comes near any other edge of the body.

## Proofs before return

Trims stay inside their edges; every new or trimmed curve on a changed flat face is proved apart from the face's
other curves (Lipschitz balls for lines, circles and ellipses; control hulls for Bezier pieces); joints neither fold
back nor cross again (closed-form line/conic roots); loops still nest (one outer loop, holes winding once); no other
body edge enters a blend's cross-section region; then `topology::check`. Work limits refuse; nothing is assumed.

## Sources (public mathematics)

- Rolling-ball constant-radius blends: the ball centre lies on both offset surfaces; for two planes the offset planes
  meet in the axis line (setback r / tan(phi/2)), for a plane and a coaxial cylinder in the tube's centre circle.
- Plane sections of a circular cylinder: circle (square plane) or ellipse c + u cos t + v sin t with u = r e1,
  v = r e2 + lambda d (the part across the axis is the circle; lambda puts it in the plane).
- Two equal circular cylinders with crossing axes meet in two plane ellipses, in their planes of mirror symmetry.
- Sine rule for the distance-angle chamfer: s' = s sin(theta) / sin(theta + phi).
- Exact volumes in the tests: Pappus's centroid theorem, the prismatoid formula and Steiner's formula for
  parallel bodies of a convex polyhedron (V + S r + (1/2) sum(edge length x exterior angle) r^2 + 4/3 pi r^3).
- Clearance proofs: Lipschitz bounds |C(t) - C(m)| <= max|C'| |t - m| and the convex-hull property of positive-weight
  (rational) Bezier curves under de Casteljau subdivision; winding numbers by summing principal angles over pieces
  whose bounding balls exclude the point.
