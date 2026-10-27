# Boolean operation boundary

Owner: Astra/Codex. Independent implementation; no OCCT or Truck code/dependency.
The public entry is `cut(body, tool, tolerance) -> Result<Body, KernelError>`.
Inputs and tolerance use the same linear unit (the worker protocol currently
declares inches). Neither input is mutated. Success returns the entire result;
failure returns no partial geometry. A completely removed body is successful
empty geometry. Unsupported geometry is never substituted, offset or repaired.

## Supported domains

### Rectilinear subtraction

Closed solids with axis-aligned planar faces and axis-aligned straight edges.
Faces may have multiple closed loops, including holes. Rectangular through cuts,
blind pockets, flush notches, exact contact, internal voids and repeated cuts are
supported. The result must contain one connected material component or be empty;
severing cuts that produce multiple solids are currently unsupported.

All input vertex coordinates define the arrangement planes. Membership is
constant inside the resulting rectangular cells for the admitted geometry. The
implementation verifies the complete operand boundaries against that membership
in all three directions. It subtracts memberships and constructs shared analytic
planar faces, edges and vertices on the surviving boundary. This is not voxel
approximation or tessellated-solid substitution. Coplanar cell seams remain.

Limits are explicit: at most 96 distinct coordinates per axis, 200,000 cells,
20,000 input faces / 80,000 edges / 80,000 vertices, and bounded loop/boundary
validation work. Distinct coordinate planes closer than the supplied tolerance
are rejected; they are not merged. Planar normals/edge directions must be exactly
cardinal in the input representation. Oblique geometry is unsupported.

### Contained circular through-bore

A plain axis-aligned solid box with six faces, twelve edges and eight vertices;
the complete box topology and geometry are verified, not inferred from its bounds
alone. The cutter is a closed axis-aligned cylinder represented by two full-circle
cap edges and one seam. It must span the box along its axis and have positive
side-wall clearance greater than the supplied tolerance.

The result preserves the box and adds circular inner cap boundaries and an
inward-facing cylindrical wall with its seam counted twice. It retains analytic
circles/cylinder rather than a polygonal approximation. Supported axes are X/Y/Z.
The original plain-box admission remains. Other closed inputs can use the
local planar passage domain below. Blind round pockets, side-wall contact/
crossing and fillet-crossing cuts remain unsupported.

### Circular passages through planar regions of curved solids

A full cylindrical tool may cross exactly one material interval bounded by two
axis-normal planar patches. Both patches must contain the entire disk, with
clearance from every outer and inner boundary greater than tolerance. Analytic
line/circle winding and minimum-distance checks establish containment. Analytic
surface-patch bounds must exclude EVERY other boundary patch from the finite
cutter corridor; inconclusive bounds reject the cut. Thus a clear web hole in
a filleted I-beam is admitted without changing its fillets. Separated repeated
bores are supported; overlapping, tangent or identical existing holes are not.

Accepted surrounding geometry consists of planar faces with line/arc loops and
simple axial cylindrical patches with two circular trims and two generators.
Checks cover finite coordinates/frames, references, closed oriented topology,
curve-on-surface consistency, planar loop orientation and cylindrical patch
orientation. This is a local cut operation on an input solid, not a general
self-intersection repair or global B-rep certification service. Work limits
remain 20,000 faces and 80,000 edges/vertices. Output clones the existing boundary
and appends two circular cap loops and an inward cylinder; unrelated geometry
is retained exactly in its original representation.

### Contained capsule through-slot

A plain axis-aligned box and a six-face capsule prism cutter: two semicylindrical
walls, two straight planar walls and two caps. The end radii must be equal, the
sweep axis cardinal, and the cutter must span the whole box thickness. The slot
may rotate within its cap plane. Its analytic bounds must have side-wall clearance
greater than tolerance. Conservative full-circle bounds may reject some slots
whose actual semicircle alone would fit; there is no guessed acceptance.

Admission matches the entire oriented cutter B-rep to a capsule prism, allowing
equivalent indexing, edge directions, loop starts and surface frames. The output
uses the ORIGINAL cap curves as a hole in the box cross-section, swept through
the box thickness by our own extrusion code. Cap curves and cylindrical walls
remain analytic. It does not replace an unsupported cutter with a capsule.
Blind, side-wall touching/crossing and repeated/interacting curved cuts remain
unsupported for slots; capsule cuts still require a plain uncut box as their base.

## Admission, errors and representation

The public entry selects a supported input domain BEFORE execution: tools with
one cylindrical face enter bore admission, two enter capsule admission, none
enter rectilinear admission, and other cylinder counts are unsupported.
An error never triggers another algorithm. A malformed model reports InvalidInput;
unsupported geometry, resolution or resource requirements report Unsupported.
The worker catches unexpected faults separately. Validation includes oriented
edge occurrences, connected loops, face orientation, and (for rectilinear solids)
surface coverage and vertex fans. Boundary curves and surfaces remain analytic,
but numerical calculations use floating-point arithmetic and declared tolerance;
this is not an exact-arithmetic claim.

The result is experimental geometry. The GCmaker adapter may inspect/display it
but must not publish it into normal saved-part or manufacturing state. Actual
native-reader compatibility and native volume checks belong to that adapter's
integration tests; Rust-only tests do not prove Eyeshot display compatibility.

## Verification

`cargo test --lib boolean` runs focused geometry and protocol checks. These
include independently calculated volume, a divergence-theorem integral of the
returned bore and slot boundaries, material probes, 60 deterministic box-cut cases,
scale/translation, existing face holes, repeated rectangular cuts, exact and
near contact, unsupported cases, malformed inputs and empty results, slot
rotations and all three sweep axes, and curve/index/loop ordering changes. The
foundation integration tests additionally exercise protocol/neutral-format
interoperability. The local source commits introducing capabilities are 9614fe8
(rectilinear), bd93e97 (contained circular bore) and 08544a5 (capsule slot).

Methods are elementary coordinate-plane arrangements, oriented boundary
construction, circle geometry and the divergence theorem. No third-party kernel
source was translated or used. OCCT comparison, where performed by the GCmaker
integration owner, is evidence alongside independent checks, not authoritative
truth or a substitute execution path.
