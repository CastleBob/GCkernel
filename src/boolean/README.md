# Boolean operation boundary

Owner: Claude (from 2026-10-29, Bob; previously Astra/Codex). Independent implementation; no OCCT or Truck code/dependency.
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
planar passage or transverse intersection domains below. Round blind pockets on
a plain box and cuts outside those explicit domains remain unsupported.

### Circular passages through planar regions of curved solids

A full cylindrical tool may cross one to 128 separated material intervals, each
bounded by two axis-normal planar patches. Every patch must contain the entire disk, with
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
and appends two circular cap loops and an inward cylinder per interval; unrelated geometry
is retained exactly in its original representation.
Entry/exit normals must alternate along the cutter, with positive separation,
and its finite depth must span every intersected wall completely. This supports
both-wall holes in rectangular tube. Capsule and rectangular through-slots use
the same admitted intervals and attach their actual extruded profiles to each
cap pair. No geometry is constructed until all intervals are admitted; an error
returns no partial body. Repeated cuts preserve the previous holes and fillets.

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

### Planar tools on curved stock

A verified plain axis-aligned box spanning exactly one complete stock coordinate
uses the clear-passage or transverse intersection domain, selected from its
footprint before construction. A clear passage retains unrelated geometry.
Box tools spanning two or more stock coordinates enter transverse cell admission.
Other planar-prism representations select the convex-profile passage below when
their finite corridor avoids cylindrical patches, otherwise transverse admission.
Their normals and the stock cylinder axes select a compatible cardinal section
chart before any operation runs. The entire cutter is still checked through its
actual oriented boundary; its shape is not inferred from its bounds. No failed
admission or operation retries another chart or solver.

This admits square tube/pipe end trims and compatible tilted cuts, including
45-degree miters produced by our extrusion. Coincident planar caps are allowed;
empty results are successful, while disconnected material is refused. The
captured native 45-degree slab currently has a tiny nonzero third component in
its cap normals and is explicitly unsupported: no automatic near-axis adjustment
is made. Exact compatible wedges are a separate admitted representation.

Checks cover web/root rectangles, repeated pipe/tube trims, independent centroid
volume integrals, 13,125 material probes per end case, transformed/scaled cases,
complete removal, flush contact and severing rejection. Square trims commute
with separated wall holes in the tested pipe and rectangular tube. General
oblique geometry remains outside this
entry's supported domain (a box spanning no stock coordinate: see "Survey gaps closed 2026-10-31"). Tests can export exact end-cut inputs/results to a
caller-selected directory with GCKERNEL_ASTRA_END_EXPORT for native review.

### Prism section cuts (Claude, 2026-10-29)

`section.rs`. Selector (`section::selects`, from operand geometry before any
construction; one domain runs and its failure is returned, never retried):
1. both operands are right prisms along the same cardinal axis k: every face is a
   cap plane normal to k (within 1e-12 of cardinal) or a plane/cylinder parallel
   to k, the caps lie on two levels, and the prism proof below holds;
2. the tool spans the body along k (tool low <= body low + tol, tool high >=
   body high - tol); the first axis k = 0, 1, 2 meeting 1-2 is the chart;
3. and ANY of: the tool's section meets or crosses the body's section boundary
   (outline or a hole); the body is no longer a plain box (six cardinal planar
   faces), e.g. Bob's p9999 hole after its corner cut; or a representative point
   of the tool's boundary (its first edge's start) is not classified inside the
   body (winding not 1), i.e. a tool wholly enclosing the body or wholly outside
   it with no boundary meetings. Only a clear interior cut in a plain box keeps
   the domains below (bore, passage, rectilinear).
Bounds: at most 2000 section edges per operand and 8000 kept pieces.
Prism proof: the section's area times the extent equals the measured volume.
The result is our own extrusion of the 2D regularized difference of the sections
(line and circular-arc loops): every edge split at closed-form meetings and at
endpoints lying on the other boundary; meeting points within the tolerance are one
vertex; coincident pieces decided once per pair by their tangents (same way:
removed; opposite: kept); every other piece kept or dropped by its midpoint's
winding; collinear lines and co-circular arcs merged so no seams remain. Refused:
a vertex met twice (a pinch, e.g. a hole tangent to the outline: non-manifold),
disconnected material, a tool not spanning the body, conic section edges.
Public method: subdivide-and-classify (Martinez, Rueda, Feito 2009), written
independently. Tests: corner triangle on the plate faces (exact, +-1 ulp), arc
notch both ways round, cut across a hole, re-cutting a hole or slot (no change),
crossing bores and slots, Bob's p9999 corner cut, all against analytic volumes.

### Radial bores through a pipe wall at any clock angle (Claude, 2026-10-29)

`radial.rs`. Selected (after the section domain) only when EVERY condition holds,
clearances included, from the geometry before anything is built: the tool is a
full round cutter (one cylinder, two caps, two full-circle cap edges); its axis is
perpendicular to a pipe's axis (outer cylinder face facing out, one coaxial inner
face facing in); hole radius plus axis offset is below the inner radius; its start
cap lies beyond the outer wall along the cutter (|(s - p).d| > Ro, a cap disk only
moves away from the axis sideways) and its end cap is either in the bore
(sqrt(((e - p).d)^2 + (offset + rho)^2) < Ri) or past the far wall ((e - p).d > Ro);
each wall is ONE face whose boundary boxes (ends, seam, earlier holes) are clear of
the new loops' box and whose end rings lie beyond the loops along the axis; every
other face's boundary boxes are clear of the cutter's box. Otherwise the earlier
domains keep the cut. Result: the body unchanged except, per pierced wall, an inner
loop on the outer and inner faces (exact CylinderPair intersections, closed, one
vertex each), one straight seam and the cutter's cylinder face between them. Pair
edge boxes: Q - c is within rho along b and w = b x a, and the point's part along a
is branch sqrt(r^2 - w^2). Tests: Bob's Pipe Drill Test holes at 0, 340 and 320
degrees from his placement (surface + 0.001, 1 deep) each remove the exact one-wall
integral; a 37-degree hole through both walls removes two; ends, seam, overlapping
holes and wall pockets are refused or left to other domains.

## Admission, errors and representation

### Transverse analytic end cuts

The captured I-beam rat-hole pair and partial-depth cope profiles now use the
same transverse cell engine described below. The earlier single-rectangle
section/subtraction routines remain only as independent test references.
Input is reconstructed from its CURRENT oriented B-rep; it does not regenerate
an uncut stock recipe. Interior openings and finite caps use complete material
cells, so previous holes and partially cut walls are retained.

Admission checks every surface against its analytic boundary, all oriented
sections and transverse interfaces. No input patch may disappear below partition
resolution. A finite end skin is never silently treated as a flush cut.
Coordinate coincidence and shared vertices use the declared B-rep tolerance;
original supporting surfaces remain analytic. Assembly splits shared boundary
edges before enforcing closed occurrence topology. Coplanar partition seams
remain in the result. Limits include 256 stations, 5,000 input faces / 20,000
edges / 20,000 vertices and explicit recognition/assembly work bounds.

Verification includes both GCmaker-captured rat-hole tools in order, type-3
conic replay, 630 production material probes, all six cardinal frames, scales
0.1 and 25.4 with translation, input preservation, malformed/unsupported cases,
and returned-boundary volume against an independent dimension-based formula.

The public entry selects a supported input domain BEFORE execution: tools with
one cylindrical face and full-circle seam topology enter bore admission;
other one-cylinder tools enter transverse cell admission. Two-cylinder tools
with distinct stored radii enter that same general profile domain without
changing either radius; equal-radius tools enter capsule admission. Zero-cylinder
tools on curved stock, and polygonal tools with oblique planar walls, enter the
planar tool domains described above; remaining planar tools enter rectilinear
admission. Other cylinder counts are unsupported.
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

`endcut_cross.rs` uses a bounded arrangement of analytic Y/Z coordinates in each
X band. Input face intervals come from oriented B-rep crossings. Winding in BOTH
section directions must agree, with each boundary present exactly once. Whole-
band graph ordering proves that the cell membership remains valid throughout
that band. Only identical graph representations share a coordinate; distinct
ones closer than the input tolerance are refused, never shifted together.
Subtraction keeps stock cells outside the finite cutter. Internal cell faces
are omitted before constructing boundary curves. Connectivity includes shared
faces within each band and positive-area overlap between bands. Existing holes
remain; disconnected output components remain unsupported. Input interface
coverage, analytic surface checks and closed oriented topology are required.

X stations come from vertices, circle extrema and analytic graph crossings;
values within tolerance are one station. Such a station is placed AT the
analytic value it agrees with: the extremum (existing), or the crossing when
all crossings within tolerance of it agree within roundoff (Claude, 2026-10-29;
Bob's rule: what agrees within the tolerance is the same). An authored vertex a
few ulps short of a crossing therefore no longer starts a band past it, which
refused a second miter on stock already mitered at its other end (Bob's HSS2x2
36 Miter Both Trapezoid and L3 Angle 48 Miter Both Ends; tests/miter_pair.rs
checks each removed volume against its exact section integral). Stations that
end up within tolerance of each other are refused. The ordering proof itself is
unchanged: a sub-tolerance reversal inside a band is still refused.

The graph ordering bound checks all candidate extrema of a difference of two
line/circle graphs. Equating their derivatives gives a closed-form critical point
for line/circle and two linear candidate equations for circle/circle. An interior
crossing or tangency without a partition event is refused; zero width at a band
endpoint is admitted only when the interior is strictly positive. A collapsed
void interval meeting an exterior X boundary is refused as nonmanifold contact.
The arrangement has at most 129 coordinates per direction, 4096 grid entries,
128 retained cells per band and two million charged arrangement operations per
cut, alongside the existing global patch/assembly limits.

Validated round/capsule tools whose finite clearance corridor meets a
perpendicular cylindrical patch enter this domain before a cut runs. Bodies
already containing cylinder-pair edges also enter it. Inconclusive classification
may reject an otherwise feasible cut; a failed operation never invokes another
solver. Supports cardinal perpendicular cylinder axes, provable rectangular
sections and positive intersection discriminants throughout each edge. Finite
cut caps can lie inside the material when the complete cell arrangement is
proved: this admits the surveyed tube-corner cut and bounded pipe/I-beam blind
pockets. Tangent branches, unresolved crossings within a band and multiple
result solids remain unsupported. At most 128 intervals per band and 4096 total
intervals are permitted, alongside the existing transverse work limits.

Cylinder/cylinder edges remain analytic `Curve::CylinderPair` (neutral type 4),
parameterized on one cylinder with the algebraic axial root on the other.
Supporting cylinders, branch and angular range survive worker replay. There
are two possible parameter charts for each intersection edge. Construction
chooses a chart whose discriminant stays positive over the entire trim, so a
large hole can cross the vertical end of a smaller fillet without a singular
parameterization. If neither chart is smooth, the cut remains unsupported.
This is an analytic coordinate choice before edge construction. There is no
polygonal geometry substitution. Their mass integrals are numerical
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
finite tool and all supporting section graphs can be ordered over each band.
Partial wall depth with provable intersection stations is admitted; unresolved near coincidences still refuse. Independent pipe
wall integrals, material probes, protocol
replay, six frames and far-translated/scaled cases verify this domain.

The tube-corner cut and two curved blind pockets are also checked against
independent piecewise dimension integrals, 10,925 material probes per case,
repeat replay, three axis placements and scales at 1e5 translation. The previous
outer-edge tangent refusal stays a regression. A 21-hole pipe sequence checks
all preceding openings via genus and accumulated volume at every protocol step.

Circular extremal stations retain their exact parameter 0 or pi when their
stored X coordinate equals cx+r or cx-r. Re-subtracting cx can otherwise yield
r minus one floating-point unit and invent a nonzero radial root. This arithmetic
fix does not snap near stations or change a supporting circle; a scaled,
far-translated miter regression exercises the original loss.


## Active graph intersection stations

Before ordering section coordinates, the engine collects active supporting
lines and circle branches from both operands. Closed-form line/line,
line/circle and circle/circle intersections add stations within each existing
vertex/extremum band. Candidate roots from squared equations are checked against
both stored branches. Interior tangencies are refused. Subsequent whole-band
ordering and all B-rep coverage checks are still required. Discovery shares the
two-million-operation budget; station and patch limits remain in force.

At a computed station, separate graph evaluations can leave endpoint residuals
of a few floating-point units. Endpoint ordering admits at most
min(input tolerance, 16*machine epsilon*sum of evaluation term magnitudes).
Interior critical points still require positive separation greater than the
input tolerance. This is a bounded endpoint arithmetic allowance, not outward-
rounded interval certification. Distinct near-coincident graphs are still
refused; no supporting surfaces are snapped or merged.

The captured Rect Tube now replays its first four tools; STR-1 replays its first
partial-depth cope. The actual tube miter has an independent dimensional
area/first-moment integral and 47,027 material probes. Later captured steps are
checked against Claude's independent Eyeshot mesh references (absolute allowance
0.001 in^3 for tessellation; see tests/fixtures/parts/*/reference.md). Native
returned-body acceptance is recorded separately by the GCmaker owner.
STR-1's next wedge retains an explicit near-coincident-line refusal. Rect Tube's
fifth, arbitrary side-profile cut uses the local passage below; the whole-body
chart still cannot represent all three cylinder axes together.
A pipe hole whose finite bottom crosses the inner wall is independently checked
by a piecewise wall integral, 10,925 material probes, repeat cuts and transformed
scale/translation tests; its volume is 70.48824992315208 in^3.


## Convex planar profile passages and exact coplanar patch union

A bounded, straight-edged convex prism can pass through planar wall intervals
without changing remote curved faces. Two exact cardinal caps and perpendicular
side planes select its axis. The complete input B-rep must match our own
extrusion of its original cap lines; the bounding disk is only a clearance
proof, never the removed shape. Three to 32 profile edges, one loop per face,
and the existing passage work limits are enforced. Examples include triangular
punches through a flat plate or a filleted web and a side cut after top holes.

Before this local passage proof, connected axis-normal cap patches on EXACTLY
the same stored plane and with the same outward direction are joined. Only
shared, opposite topological edge occurrences are cancelled. Surviving edges
must assemble into one outer loop and resolved holes without branching. Unused
seam edges/vertices are removed and references remapped; analytic surfaces and
surviving curves are unchanged. Near planes are never merged. Oriented analytic
validation runs before and after the union. This is a deterministic exact-plane
boundary union, not a retry after an unsuccessful Boolean or tolerance healing.
The entire finite passage is admitted before cavity construction; every actual
profile line is extruded through the proven wall depths and spliced into the
joined caps. Failures return no partial body and never mutate the input.

All FIVE captured Rect Tube tools now replay, final volume54.34233932006709in^3.
The fifth cut removes the independently computed triangle area times the0.25in
wall depth;15,283 independent local ray/material probes agree. Union preserves
volume/area, is idempotent, keeps near planes separate and rejects malformed
input. Previous non-cap faces remain unchanged after reference remapping.
The fifth cut is tested on three axes and at1e5 translation. General convex
profile checks cover a plain plate, two separated I-beam web triangles on all
three axes, repeated-hole and partial-depth refusal. Native completion remains
a separate handoff result, not implied by these Rust tests.


## Survey gaps closed 2026-10-31 (Claude)

Bob's linear-nesting parts that the GCmaker native survey found stopping in the kernel. Fixtures are the worker's
own inputs (tests/fixtures/survey_gaps); tests/survey_gaps.rs checks every removed volume against an independent
value from the parts' dimensions (closed form, or a quadrature of the analytic section written in the test).

- **Box tools spanning no stock axis** (L3 Angle copes): a plain box that spans none of the curved stock's axes
  now enters the transverse chart (planar_cross_cut) instead of being refused. A blind web pocket on the I-beam
  takes the same path (exact 0.92 in^3 in box_cut's test, formerly a refusal). The chart is chosen from the tool's
  normals AND the stock's planes: a sweep w whose beam axis v and w would both be in an earlier miter's normal has
  no axial span for that plane, so the chart is not chosen (L3 Angle 18 Miter And Cope; geometry before the cut).
- **Coplanar cell seams**: an earlier transverse result keeps partition seams; section segments on the IDENTICAL
  graph, facing the same way and meeting at the seam's crossing, are one boundary segment and are joined (L3 Angle
  8 Coped Both Ends). Distinct graphs, opposite orientations or gaps are left alone.
- **Analytic extremal stations**: a circle whose plane holds the station direction (a cylinder across x) gives its
  extremum station as c.x +- r, not its evaluation through its own frame. A wall arc left by an earlier cut,
  evaluated at its extremum, gave 0.8299999999999998 for cx + r = 0.83; acos(1 - e) ~ sqrt(2 e) turned those two
  ulps into 2e-8 of angle, an edge end 1e-8 off its vertex and 7.6e-8 in^3 of volume (Pipe 1.25 18 Copes Rotated
  90 and the slot after a miter, both now exact to 1e-10 against their integrals). A snap of circle_angle itself to
  the extremum within roundoff was tried and rejected: at a genuine crossing a few ulps from an extremum it moved
  the rat-hole tests' pair ends and lost 3.7e-7 against their independent references.
- **Radial bore on the wall's seam** (Pipe 1.25 12 Cross Holes, a side hole at 3:00): the seam of a plain periodic
  wall (one full ring at each end, one straight seam used both ways, nothing else at its vertices) is a
  parametrization choice, not geometry; when the new hole's angular interval reaches it, the seam moves to the
  middle of the widest gap between the wall's holes (rings restarted at the new generator, surfaces unchanged).
  Each hole's interval: centre direction branch * a, half-width asin((|offset| + rho) / r). The radial domain
  now also refuses planar faces turned inside out (outward signed area), which the topology check alone misses.
- **Split transverse cuts** (endcut_split.rs; Pipe 1.25 18 Copes Rotated 90, Pipe 1.25 24 Slot And Miter): when the
  body has faces the chart cannot hold (a cylinder along the stations, a plane turned in both v and w) that lie
  wholly beyond the tool along the beam axis, the body is split by y = c in the middle of the widest gap between
  vertices and edge extrema there; every face y = c crosses must be a plane with n.y = 0 or a cylinder along y,
  crossed by exactly two straight edge uses. Each piece gets a cap at y = c; the near piece is cut in the chart;
  the caps are removed and the open edges along y = c are split at each other's vertices and joined where they are
  the same curve (one vertex within tolerance), every one of them or the cut is refused. With two longitudinal
  cylinder axes besides the sweep, the beam axis is the one for which the other axis's faces are separable this
  way (geometry before the cut). The far piece keeps its faces, curves and surfaces exactly; the faces y = c crossed
  keep a seam there.
- **Oblique pipe copes** (cope.rs; Pipe 1.25 12 Coped 60 To 2.375): a full round cutter oblique to a plain pipe's
  axis that holds one whole end annulus strictly inside its cylinder, clears the other end's whole disk, and whose
  finite length covers every shared point on the material side. Each wall's new end is one CylinderPair loop (wall
  = parameter cylinder, cutter = root, branch toward the material), discriminant proved positive all the way round;
  the cutter's cylinder between the loops faces into the cope. Bounds: with phi between the axes, a shared point is
  O_b + alpha b + beta d + gamma n with |alpha| <= Rc / sin phi, |beta| <= Ro / sin phi, clipped by the coped
  end's half-space; its place along the cutter is extreme at that clipped box's vertices. The squared distance from
  the cutter's axis round a circle across the pipe is a degree-two trigonometric polynomial (eight samples, chord
  bound). Anything else (a pipe already holed or cut, a short cutter, an end not wholly inside) is refused by name.

Not done in that round: native normalized ends (closed below) and NC1 rolled-root cutters.

## Normalized pipe ends (Claude, 2026-10-31)

Bob's "short pipe test" (7 x 0.5 pipe, Normalize End over a 13 diameter Preserve Length cope). GCmaker's rule
(PipeCopeNormalizer, EndConditionSideCompositeBuilder): a radial torch cuts one boundary per angle, so at each angle
the end keeps the deeper of the two walls' crossings with the mating cylinder and carries it radially across the
wall. The exact face is a right conoid with radial rulings, `Surface::Conoid` (wire surface 9):
S(t, s) = c + s rho(t) + h(t) a, h the deeper wall's CylinderPair height; its rims are `Curve::RadialPair` (wire
curve 7: a CylinderPair's height carried to radius s). Construction (`normalize.rs`, worker op 19) proves the deeper
wall: both walls' discriminants positive all the way round (else "does not cross both pipe walls all the way round")
and h_in - h_out of one sign within the tolerance at every angle, by samples less max|D''| step^2 / 8 with |h''|
bounded from the pair's trigonometric discriminant (`geom::pair_bend_bound`); a square cope's two crossings meet at
the pipe's sides, which is the same within the tolerance. A wall that changes round the pipe (tilted or off-centre
copes) is refused by name ("the deeper wall changes round the pipe"). The tool is GCmaker's waste band (radii sa..sb
about the pipe, far cap beyond the end, conoid top). `normal_end.rs` cuts it from a plain pipe on the same axis line
when the band reaches past both walls, the far cap lies beyond one end, and the conoid's height is proved strictly
between the ends (else "the envelope reaches a pipe end"); the far cap, walls and far rings are kept exactly, each
wall trimmed to its rim. Bodies with a normalized end take only the radial bore domain, whose clearance bounds a rim
in 128 pieces (circle-piece extrema plus the proved height range), the pieces at the hole's own angles giving that
end's place along the axis; anything else is refused by name. Mass: X.N = -s (C.a + h) + h' C.v with Green's theorem
in (t, s); area from the closed-form antiderivative of sqrt(s^2 + h'^2) in s (quadrature in t). Tests
(tests/normal_end.rs): Bob's end against (Ro^2 - Ri^2) / 2 times the integral of y0 + sqrt(Rc^2 - Ri^2 sin^2 t)
(1e-9), its face area, the tool's own volume, his top hole after it (one-wall integral, 1e-10), the far end, a turned
cope, the display mesh within area x chord, op 19 and the neutral round trip, and the named refusals.
