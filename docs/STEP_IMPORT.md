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

Assemblies/occurrence transforms/mapped items, B-spline curves/surfaces, pcurve
reconciliation, arbitrary swept/revolved supports, void shells, faceted/non-solid
models, unsupported shape representations and multiple/named DATA sections are
refused explicitly. Old code-page shifting string directives are not supported.
IGES/STL and STEP export are not part of this increment. Metadata beyond body
names (colors, PMI and product-tree semantics) is not returned.

Limits: 64 MiB input, 500,000 entities, 4,000,000 values, nesting 64, unit-chain
depth 32, string storage length 32767 excluding quotes, positive source IDs no
larger than 2^53-1. Errors reference entity IDs or filtered-stream byte offsets.

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
