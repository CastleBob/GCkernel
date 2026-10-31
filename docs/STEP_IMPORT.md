# Independent STEP intake

Owner: `step` (Codex, Bob's 2026-10-31 instruction). No third-party dependencies
or implementation source. Public references: ISO 10303-21 grammar at
https://steptools.com/stds/step/IS_final_p21e3.html and EXPRESS entity descriptions
at https://www.steptools.com/docs/stp_aim/ (particularly edge_curve, face_bound,
axis2_placement_3d, surface_of_linear_extrusion and surface_of_revolution).

## Boundary

`step::read` / `read_file` returns all admitted `ImportedBody` values or one error.
Input is never written. The module owns Part 21 parsing, representation-scoped
units, entity resolution and exact conversion to the existing neutral Body.
Geometry/topology/mass/tessellation/boolean algorithms retain their owners.
No repair, sewing, endpoint snapping, curve fitting, alternative importer or
partial-success result. Source IDs are file-local, not application reference names.

## Initial supported domain

- Part 21 HEADER plus one unnamed DATA section, simple/complex entities, forward
  references, comments, enumerations, typed values, quoted/escaped strings,
  X/X2/X4 and UTF-8. Physical control octets are ignored per section 5.2, including
  wrapped keywords/strings. Duplicate/missing references and truncation fail.
- AP203 CONFIG_CONTROL_DESIGN, AP214 AUTOMOTIVE_DESIGN and AP242 MIM schema names
  are recognized. This is an entity subset, not complete AP/schema validation.
- Explicit MANIFOLD_SOLID_BREP / CLOSED_SHELL / ADVANCED_FACE / EDGE_LOOP graphs
  in standalone SHAPE_REPRESENTATION or ADVANCED_BREP_SHAPE_REPRESENTATION.
  Each face requires one explicit FACE_OUTER_BOUND; inner FACE_BOUND loops remain
  exact. Multiple standalone solids are retained separately, in source-ID order.
- LINE, CIRCLE and ELLIPSE edge supports, trimmed by EDGE_CURVE vertices and
  same_sense; oriented edges/bounds and face sense are applied independently.
- Clamped polynomial or positive-weight rational B_SPLINE_CURVE_WITH_KNOTS of degree 1..3, including
  complex entity encoding, over the full domain or explicit parameter trim below. Exact
  homogeneous knot insertion and degree elevation produce cubic Bezier spans; no fitting.
  Both incident faces share each decomposed edge and intermediate vertex.
  same_sense and oriented-edge reversal reverse the complete span chain.
- Plane/cylinder/cone/sphere/torus supports when their boundaries are admissible
  to existing GCK validation. Exact line extrusion -> plane, normal circle/ellipse
  extrusion -> cylinder/elliptic cylinder, meridian-line revolution -> plane,
  cylinder or cone. These are analytic representation changes, not approximations.
- Per-representation SI and positive conversion-based length/angle units. Output
  inches, radians, cubic inches. Directions and axis-placement derived frames
  follow the format definition. Unit magnitude is never inferred from the part.

## Admission and limits

Strict neutral decoding, connected closed-shell topology, curve/support agreement,
planar outer/inner orientation, positive per-curved-face area and positive finite
total volume/area are required. These checks do not establish a general geometric
self-intersection detector or healing system. Import success does not guarantee
all later boolean domains support the body. Tessellation remains display only.

Assemblies/occurrence transforms/mapped items, higher-degree/unclamped
splines, Cartesian-only/preferred or nested spline trims, general spline surfaces, pcurve
reconciliation, arbitrary swept/revolved supports, void shells, faceted/non-solid
models, unsupported shape representations and multiple/named DATA sections are
refused explicitly. Old code-page shifting string directives are not supported.
IGES/STL and STEP export are not part of this increment. Metadata beyond body
names (colors, PMI and product-tree semantics) is not returned.

Limits: 64 MiB input, 500,000 entities, 4,000,000 values, nesting 64, unit-chain
depth 32, string storage length 32767 excluding quotes, positive source IDs no
larger than 2^53-1. Errors reference entity IDs or filtered-stream byte offsets.
Spline curves are limited to 8192 controls and bodies to 100000 decomposed
edges. Bernstein surface bounds allow 8192 subdivision visits and depth 24;
failure to establish the bound is a refusal, not sampled acceptance.

## Polynomial boundary admission (2026-10-31)

Planes retain their control-hull bound. Circular/elliptic cylinders additionally
admit Bezier boundaries using degree-six Bernstein bounds on the squared radial
residual, tightened by subdivision with a floating-point roundoff allowance.
This checks the complete span, not a selection of sample points. It is not a
general interval-arithmetic certification or arbitrary self-intersection test.
Cylinder mass uses the existing analytic surface chart and boundary quadrature;
MassProperties.exact is false for that calculation. No source curve is projected,
snapped, approximated or rewritten to agree with its surface. Other Bezier face
domains retain their existing acceptance/refusal policy.

The worker/bridge default remains 1e-6 inches. File-declared uncertainty is not
silently substituted for that limit. The seven previous spline corpus failures
now pass decomposition and fail geometric admission at this tolerance. Examples
of sampled departures on the first refused span: Ellipse Pipe 8.90e-6 inches,
TubePipe 1.02e-4, Pipe_10in_Slot_3holes 5.06e-5. These are sampled witnesses, not
whole-file maxima. The assembly also contains spline surfaces and stays refused.
The 19-file default admission count therefore remains 11, plus a separately
authored polynomial-cylinder fixture that verifies the new supported domain.

Developer-only comparison: `step_inspect --tolerance <inches> <files...>` makes
one explicitly requested read per file; it never retries or changes app policy.
An experiment at 0.0004 inches admitted Ellipse Pipe and Pipe_10in_Slot_3holes;
TubePipe still failed (another span departs by about 0.00183 inches). This is
not approval to raise the production tolerance or a claim of manufacturing parity.

Mathematical references: C.-K. Shene's MTU notes on knot insertion/subdivision,
https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/multiple-time.html
and https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/subdivision.html;
Bernstein basis properties, https://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node9.html.
The independent Cox-de Boor test oracle compares every span for degrees 1..3 and
all supported internal multiplicities. Further tests cover shared topology,
reverse sense, units, analytic mass bounds, malformed data, a hidden radial bulge,
neutral roundtrip and tessellation. GCm's native comparison uses an independently
authored fixture, with source preservation; no CAD library becomes a runtime dependency.

## Rational boundary milestone (2026-10-31)

RATIONAL_B_SPLINE_CURVE weights must be finite, positive, one per control point,
with largest/smallest ratio at most 1e12 (explicit numerical-domain limit).
Weights are scaled together before homogeneous knot insertion/degree elevation;
this leaves the rational function unchanged, including its parameterization.
The resulting cubic has neutral curve tag 6: `[6, p0.xyz, p1.xyz, p2.xyz, p3.xyz,
w0, w1, w2, w3]`, 17 doubles total. This is distinct from surface tag 6.
The op-16 envelope and default tolerance do not change. Standalone transforms
and body roundtrips retain the emitted weights. Old workers cannot consume tag 6.

Rational boundaries are admitted on planes and circular/elliptic cylinders.
Positive weights give the planar control-hull bound. For cylinders the degree-six
Bernstein numerator `Q.Q - radius_squared * W*W`, divided by a positive lower
bound on `W*W`, bounds the whole radial residual. Roundoff/work-limit refusal is
explicit; this is not formal interval arithmetic. Mass uses boundary quadrature
and reports `exact = false`. Display subdivision bounds the positive-weight
control hull's distance from the chord; depth 24 and 200000 visits are hard limits.
Booleans and profile extrusion explicitly refuse rational boundaries. General
rational surfaces, higher degrees and unclamped/periodic knots remain outside
this milestone; a closed clamped circle is supported. Explicit parameter trims
are supported by the subsequent increment below.

Independent fixture `tests/fixtures/rational_cylinder.step` represents exact
radius-2, height-4 cylinder rings with rational quadratic quarter circles, then
exact homogeneous elevation. Tests require volume `16*pi` and area `24*pi`
within 1e-9, before/after a skew rotation and distant translation; shared topology,
neutral roundtrip, analytic derivative, weight scaling, subdivision, malformed
weights, surface mismatch, cut refusal and tessellation are also exercised.
Weighted Cox-de Boor basis tests cover degrees 1..3 and each internal multiplicity.
GCm stores the identical fixture for independent EyeShot intake comparison.

Public mathematics/specification (no implementation source):
https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/NURBS/NURBS-knot-insert.html
and https://www.steptools.com/stds/stp_aim/html/t_rational_b_spline_curve.html.

## Explicit spline parameter trims (2026-10-31)

Simple TRIMMED_CURVE entities can directly reference supported polynomial/rational
B_SPLINE_CURVE_WITH_KNOTS. Each trim set must contain one PARAMETER_VALUE, with
at most one additional CARTESIAN_POINT. PARAMETER preference selects the given
parameter; a redundant point must agree within the unchanged request tolerance.
UNSPECIFIED is admitted only for parameter-only sets. Cartesian-only/preferred,
ambiguous dual-selector and nested trims are explicit refusals, never inversion
or substitution. Other trimmed basis kinds remain outside this increment.

Trim parameters are in the source knot domain, not normalized [0,1] and not
scaled by length units. Span extraction uses exact de Casteljau subdivision
(homogeneous for rational curves); all interior knot spans retain shared topology.
SenseAgreement directs the trim from trim_1 to trim_2; EDGE_CURVE.same_sense and
ORIENTED_EDGE/FACE_BOUND senses remain separate. A closed clamped basis can be
traversed across its seam only when declared closed and its seam agrees within
tolerance. This is not support for unclamped periodic knot vectors. Out-of-domain,
equal/zero-length or direction-inconsistent open trims fail explicitly. Trimmed
endpoints must agree with both EDGE_CURVE vertices; the source curve is never
snapped or fitted to those vertices. Bare-spline vertex-only interior trimming
still requires a future unambiguous point-inversion implementation.

Evidence: independent weighted Cox-de Boor oracle across degrees 1..3 and all
internal multiplicities, nonuniform source knot domains, four orientation
combinations, polynomial boxes in inches/mm, malformed/contradictory selectors,
cyclic/nested basis refusal and a complete rational cylinder whose eight source
ring edges each cross a knot or the closed seam. That fixture produces 16 vertices,
17 shared edges and 3 faces, with analytic volume/area and display-mesh checks;
GCm carries an identical fixture for the independent EyeShot comparison.

Independent native comparison limitation: EyeShot agrees on the parameter-trimmed
polynomial box, rational interior trims without seam crossing, and rational seam
trims at knot boundaries. It disagrees on the combined interior/seam fixture:
Rust/analytic volume 50.2654824574 cubic inches, native mesh 5.8306674805. Uniform
knots and forward senses alone do not remove the discrepancy; an equivalent
interior trim on a shifted basis without seam crossing agrees. The combined case
is retained for Rust analytic/basis/topology/display regression, not claimed as
native agreement. This does not by itself establish the cause inside EyeShot.

Specification: STEP Tools trimmed_curve and trimming_preference EXPRESS pages;
ISO 10303-42:2003 section 4.4.40 gives the directed inherited parameter definition.
Mathematics: MTU B-spline subdivision notes cited above. No kernel source or new
dependency. Geometry/wire/bridge owners are unchanged; only STEP conversion expands.

## Vertex trims and tensor-product surfaces (2026-10-31)

Bob explicitly requested both increments. Bare spline EDGE_CURVE endpoints and
Cartesian TRIMMED_CURVE selectors now use bounded cubic point inversion. The
largest control-net coordinate defines a cubic homogeneous plane equation;
derivative roots partition it into monotone intervals, all roots are checked in
3D, and a missing or ambiguous answer fails. This is conservative inversion,
not nearest-point recovery: an ill-conditioned coordinate or an off-curve point
can be refused even when a nearby parameter exists. Closed seam endpoint aliases
are identified; the declared full-domain case retains its existing semantics.
Explicit selector preference is honored; redundant selectors must agree. Exact
subdivision, independent trim/edge/face senses and the original knot domain remain
unchanged. Nested trims, unclamped/periodic knots and guessed recovery stay out.

B_SPLINE_SURFACE_WITH_KNOTS now supports clamped degree 1..3 in each direction,
polynomial or positive-weight rational, through independent homogeneous knot
insertion/elevation in U and V. The exact result is a tensor cubic patch grid,
not a mesh or analytic-surface fit. Maximum source grid 16384 controls, maximum
4096 patches, weight ratio at most 1e12. Declared closed/self-intersecting surfaces
are refused. Internal patch boundaries must agree within 1e-9 inches. Geometry
uses normalized parameter knots; STEP length units apply only to control points.

Each admitted face must cover its complete natural rectangular domain with one
loop. Lines and cubic spline edges may span multiple natural knot intervals;
endpoint interval correspondence and whole-span Bernstein difference/control-hull
bounds certify the boundary. Nonlinear reparameterizations whose equivalence is
not established by these bounds are refused. Arbitrary UV cutouts/holes, pcurve
reconciliation and periodic surfaces are not implemented. Interval bounds on the
rational derivative numerators certify a nonzero Jacobian through subdivision
(depth 12, 16384 visits per patch). This proves local regularity, not global
self-intersection freedom. Floating-point arithmetic uses outward-rounded bounds
and explicit roundoff allowances; no formal exact-arithmetic claim is made.

Neutral surface tag 8 is `[8, nu, nv, U(nu+1), V(nv+1), patches...]`; each patch
has 16 xyz controls in U-major order, then 16 positive weights. Rigid transforms
map controls and preserve weights/knots. The bridge walks this bounded record
inside bodies, with no C# spline evaluator or standalone typed surface API.
Mass integrates the entire admitted parameter domain by bounded quadrature and
reports `exact=false`. Inner integration failures and evaluation limits propagate.
Display retains canonical shared body-edge points, refines only internal edges
on both adjacent triangles, and bounds surface-minus-triangle error by degree-(4,4)
Bernstein coefficients over each triangle's parameter bounding rectangle. A
50000-triangle budget or an unproved bound is an error. Booleans and profile
reconstruction explicitly retain their existing domains and refuse spline faces.

Regression cases include independent tensor Cox-de Boor evaluation/derivatives
for all nine degree pairs, weighted restriction, ambiguous inversion, singular
charts, full solid polynomial/rational/curved-boundary/multispan fixtures,
reversed charts, transforms, shared mesh boundaries and neutral transport.
GCm's native EyeShot comparison uses the same independently authored fixtures.
The original 19-file corpus remains 11 admitted/8 refused at 1e-6 inches.
This does not establish arbitrary CAD parity or switch production STEP intake.

Public specification/mathematics, not implementation source:
- https://www.steptools.com/stds/stp_aim/html/t_b_spline_surface_with_knots.html
- https://www.steptools.com/stds/stp_aim/html/t_rational_b_spline_surface.html
- https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/surface/bspline-construct.html

## Operation 16

GCKQ v1 is unchanged. Request numeric payload: `[1, n, UTF-8 path bytes...]`,
1..32768 nonzero path bytes. Reply: `[1, bodyCount, (sourceId, nameByteCount,
UTF-8 name bytes..., sourceLengthToInches, volume, bodyValueCount,
neutralBody...)...]`. One unsupported body fails the whole operation. Ops 1-15
retain their semantics. The C# bridge only transports and checks structure.

## Verification and developer use

`cargo test --release --test step_import` runs independent analytic fixtures,
parser rejection/mutation cases, neutral round trips and a subsequent through-hole.
Set `GC_STEP_CORPUS` to GCm's `GCmaker/TestParts/STEP FILES` to exercise all 19
real samples (the corpus test returns without running the external corpus when
that variable is absent; the recorded acceptance run sets it explicitly).
`cargo run --release --example step_inspect -- path.step ...` reports bodies or
specific refusals without writing files. Exit 1 means at least one refusal.

Initial corpus: 11 admitted, 8 explicitly refused (seven spline-boundary files
and one assembly). All 11 tessellate and pass independent EyeShot body-count and
mesh-volume comparison through GCm tests (mesh volume allowance 0.05% with an
absolute floor). Several old EyeShot fixtures declare metres for inch-looking
coordinates; declared units are respected by both intake paths, never guessed.
This does not switch GCmaker's existing UI importer or establish release readiness.

## Polygonal trims on spline graph charts (2026-10-31)

An independently written affine inverse-chart certificate now admits polygonal UV
outer loops and holes on qualifying tensor spline surfaces. Homogeneous degree-(4,4)
identities prove the chart over every patch. Each line/cubic physical edge must
project monotonically onto a straight UV segment. Knot crossings split the exact
curve; degree-21 Bernstein composition residuals bound the entire edge against
its surface at the unchanged request tolerance. Inconclusive charts use their
separately defined natural-domain admission; failed trim checks never retry or heal.
Touching/intersecting/nested/outside holes and inconsistent winding are refused.
There is a 2048-edge trim budget. Curved UV trims and general pcurves remain open.

Mass uses Green's theorem over the admitted polygon, with bounded nested
quadrature, explicit knot splitting and propagated failure. Display triangulates
outer/inner loops and applies the existing whole-triangle certificate while
retaining shared body-edge points. Natural rectangular display retains its
conforming fan construction. No new neutral tag or C# geometry is required.
Independent analytic solids exercise a quadratic roof with/without a hole,
multiple knot spans, rigid frames, exact volume/area formulas, shared mesh
boundaries, transport and rejection of mismatched/reversed/duplicated loops.

## STEP serialization (2026-10-31)

`step::write` serializes a complete batch of checked bodies to deterministic AP203
text in inches. It preserves shared vertices/edges, face senses, loop orientation,
line/circle/ellipse and polynomial/rational cubic boundaries. Plane, circular and
elliptic cylinder, cone, sphere, torus and tensor spline supports are retained.
Conjugate ellipse diameters use the exact principal-axis rotation identity.
Tensor grids use degree-three piecewise Bezier knots; adjacent homogeneous control
nets must agree within floating-point roundoff after consistent weight scaling.
Cylinder-intersection curves and general revolved spline supports explicitly refuse.
There is no tessellation, curve fitting, object omission or file writing.

The strict reader validates emitted text before it is returned; this is an output
check, not an independent oracle. Names use STEP escaping (including Unicode).
Limits: 1024 bodies, 32768-byte names, 100000 vertices/edges/faces per body,
16384 spline controls, 64 MiB text and one million emitted entities.
Operation 17 request: [1, count, (nameByteCount, UTF8 bytes, bodyValueCount, body)...].
Reply: [1, byteCount, ASCII STEP bytes]. KernelStep.Write transports detached bytes;
callers explicitly own any destination write. Product/document intake is unchanged.
Verification includes analytic mass/topology, Unicode, multiple bodies, malformed
packets, strict round trips and independent EyeShot comparisons of four trimmed
roof variants and four exported spline fixtures (29 native cases total).
