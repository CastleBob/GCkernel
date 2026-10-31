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
- Clamped non-rational B_SPLINE_CURVE_WITH_KNOTS of degree 1..3, including
  complex entity encoding, whose edge uses the complete curve domain. Exact
  knot insertion and degree elevation produce cubic Bezier spans; no fitting.
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

Assemblies/occurrence transforms/mapped items, rational/higher-degree/unclamped
splines, interior spline trims, general spline surfaces, pcurve
reconciliation, arbitrary swept/revolved supports, void shells, faceted/non-solid
models, unsupported shape representations and multiple/named DATA sections are
refused explicitly. Old code-page shifting string directives are not supported.
IGES/STL and STEP export are not part of this increment. Metadata beyond body
names (colors, PMI and product-tree semantics) is not returned.

Limits: 64 MiB input, 500,000 entities, 4,000,000 values, nesting 64, unit-chain
depth 32, string storage length 32767 excluding quotes, positive source IDs no
larger than 2^53-1. Errors reference entity IDs or filtered-stream byte offsets.
Polynomial curves are limited to 8192 controls and bodies to 100000 decomposed
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
