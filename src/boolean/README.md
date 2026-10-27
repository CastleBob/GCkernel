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
The original plain-box admission remains. Other closed inputs use the local
planar passage or transverse intersection domains below. Blind round pockets
and cuts outside those explicit domains remain unsupported.

### Circular passages through planar regions of curved solids

A full cylindrical tool may cross exactly one material interval bounded by two
axis-normal planar patches. Both patches must contain the entire disk, with
clearance from every outer and inner boundary greater than tolerance. Analytic
line/circle winding and minimum-distance checks establish containment. Analytic
surface-patch bounds must exclude EVERY other boundary patch from the finite
cutter corridor; inconclusive bounds reject the cut. Thus a clear web hole in
a filleted I-beam is admitted without changing its fillets. Separated repeated
bores are supported; overlapping, tangent or identical existing holes are not.

Accepted surrounding geometry consists of planar faces with line/arc/conic loops
and simple cylindrical strips with two conic trims and axial generators. The
trims must cover the same angular interval, have nonoverlapping axial ranges,
and agree with their supporting cylinder over the entire curve. Thus separated
web holes and slots can follow the rat-hole cuts. Conic boundary-distance checks
use conservative analytic boxes; inconclusive clearance is rejected.
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
Other supported passage bodies, including the real filleted I-beam, use the same
local planar-passage proof with a disk enclosing the entire capsule. That disk
is a conservative clearance bound only; the actual original slot curves define
the cavity. Its inverted side faces and cap loops are spliced into the existing
body. This supports separated slot/hole sequences and retains existing fillets.
The enclosing disk can reject a slot that would fit in a narrower region.
Capsules meeting perpendicular cylindrical patches enter the transverse domain
below before execution; clear passages retain this local splice path. Blind and
other unsupported interacting cuts are rejected.

### Rectangular through-tools on curved stock

The cutter must be a verified plain axis-aligned box. Exactly one coordinate
axis must span the complete stock bounds; a partial-depth or ambiguous sweep
is rejected. The footprint's enclosing disk selects clear-passage versus
perpendicular-cylinder intersection admission before execution. A clear passage
splices the exact rectangular cavity and preserves unrelated geometry; a root
crossing uses the transverse plane/cylinder construction. No unsupported cut is
retried. Rotated rectangular tools and locally-through tools that do not span
the overall stock bounds are outside this bounded admission.

Tests cover repeated web slots, a root-crossing rectangle with independent
analytic volume, all three coordinate sweep axes, a rectangular slot after a
type-4 fillet circle through the wire protocol, and blind/ambiguous rejection.

## Admission, errors and representation

### Transverse analytic end cuts

A separate bounded domain handles the captured I-beam rat-hole pair, without
hardcoded I-beam dimensions. In a right-handed cardinal frame, each admitted
transverse band contains one rectangular material interval bounded by plane or
cylinder graphs. The input is reconstructed from its CURRENT oriented B-rep;
it does not regenerate an uncut stock recipe. Whole-depth cutters can remove
one end interval. The second rat-hole therefore preserves the first cut.

Admission checks every surface against its analytic boundary, the rectangular
sections between all vertex/extremum stations, and the complete oriented
coverage of transverse faces. Interval bounds prove cutter coverage across the
whole band. No input patch may disappear below partition resolution. A finite
end skin is never silently treated as a flush cut. Coordinate coincidence and
shared vertices use the declared B-rep tolerance; original supporting surfaces
remain analytic. Assembly splits shared boundary edges before enforcing closed
occurrence topology. Coplanar partition seams remain in the result.

Plane/cylinder trims use analytic circular or conic arcs, including genuinely
oblique planes. Cylinder/cylinder trims, interior transverse pockets/holes,
multiple material components, and cases without a provable rectangular section
remain Unsupported in this domain. Circular web holes still use passage
admission. End cuts currently require a cutter with one cylindrical side wall
and a base with a perpendicular longitudinal cylinder axis. Limits: 256 stations,
5,000 faces / 20,000 edges / 20,000 vertices, and explicit recognition/assembly
work bounds. This is not a general surface-intersection Boolean solver.

Verification includes both GCmaker-captured rat-hole tools in order, type-3
conic replay, 630 production material probes, all six cardinal frames, scales
0.1 and 25.4 with translation, input preservation, malformed/unsupported cases,
and returned-boundary volume against an independent dimension-based formula.

The public entry selects a supported input domain BEFORE execution: tools with
one cylindrical face and the full-circle seam topology enter bore admission;
other one-cylinder tools enter transverse-end-cut admission; two enter capsule
admission. Zero-cylinder tools on curved stock enter rectangular-through
admission; otherwise they enter rectilinear admission. Other cylinder counts
are unsupported.
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

## Transverse cylinder intersections

`endcut_cross.rs` extends the transverse construction to multiple separated Y
material intervals per X band. Input sections come from the current oriented
B-rep. Whole-depth cutters split/remove intervals; positive-area connectivity
between adjacent bands proves that one material component remains. Existing
material holes are retained. Input interface coverage, analytic surface checks,
bounded work and final oriented closed topology are still required.

Validated round/capsule tools whose finite clearance corridor meets a
perpendicular cylindrical patch enter this domain before a cut runs. Bodies
already containing cylinder-pair edges also enter it. Inconclusive classification
may reject an otherwise feasible cut; a failed operation never invokes another
solver. Supports cardinal perpendicular cylinder axes, provable rectangular
sections and positive intersection discriminants throughout each edge. Tangent
branches, blind curved cuts, unresolved crossings within a band and multiple
result solids remain unsupported. At most 128 intervals per band and 4096 total
intervals are permitted, alongside the existing transverse work limits.

Cylinder/cylinder edges remain analytic `Curve::CylinderPair` (neutral type 4),
parameterized on one cylinder with the algebraic axial root on the other.
Supporting cylinders, branch and angular range survive worker replay. There
are two possible parameter charts for each intersection edge. Construction
chooses a chart whose discriminant stays positive over the entire trim, so a
large hole can cross the vertical end of a smaller fillet without a singular
parameterization. If neither chart is smooth, the cut remains unsupported.
This is an analytic coordinate choice before edge construction.
There
is no polygonal geometry substitution. Their mass integrals are numerical
quadrature, explicitly distinguished from the closed-form plane/conic cases.
The GCmaker display adapter separately owns tolerance-bounded spline conversion.

Checks include the actual circle-over-fillet capture, all 315 production probes,
an independent dimension-based volume integral, repeated separated circles via
the public protocol, six cardinal frames, and a capsule crossing the root fillet.
Scaled/translated cases (including offsets of 1e5), nonconcentric circles, and
rat holes followed by a fillet circle, web slot and bore also pass. Native
acceptance of the complete new cuts remains separate in the shared handoff.

Annular caps and separated tube walls use oriented boundary crossings to
identify each face interval. Material sections may be separated in Y or Z;
concentric circular graphs have analytic ordering bounds across each X band.
This admits one-wall and both-wall pipe holes and pipe capsule slots whose
finite tool spans each intersected wall completely. Partial-depth wall cuts
remain unsupported. Independent pipe wall integrals, material probes, protocol
replay, six frames and far-translated/scaled cases verify this domain.
