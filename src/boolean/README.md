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
refused; no supporting surfaces are snapped or merged. One surface stored twice
is not two (Claude 2026-10-31, `same_graph`): two line graphs agreeing within 64
machine epsilons of their terms at both band ends (their difference is linear),
or two circles equal to that roundoff, are one coordinate, the first kept. This
is evaluation roundoff of a single plane (aa5's preparation wedge on the
adjoining leg's miter plane, STR-1's wedges: 1 to 2 ulps), never the tolerance:
5e-7 or 1e-9 apart is still refused.

The captured Rect Tube now replays its first four tools; STR-1 replays its first
partial-depth cope. The actual tube miter has an independent dimensional
area/first-moment integral and 47,027 material probes. Later captured steps are
checked against Claude's independent Eyeshot mesh references (absolute allowance
0.001 in^3 for tessellation; see tests/fixtures/parts/*/reference.md). Native
returned-body acceptance is recorded separately by the GCmaker owner.
STR-1's wedges, once refused as two stored lines, are one plane each (see
`same_graph`): all six captured steps now replay within 0.001 in^3 of the Eyeshot
references (final 98.79992 vs 98.80041). Rect Tube's
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
wall. The exact face is a right conoid with radial rulings, `Surface::Conoid` (wire surface 10; 8 and 9 are spline patches):
S(t, s) = c + s rho(t) + h(t) a, h the deeper wall's CylinderPair height; its rims are `Curve::RadialPair` (wire
curve 7: a CylinderPair's height carried to radius s). Construction (`normalize.rs`, worker op 23) proves the deeper
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
cope, the display mesh within area x chord, op 23 and the neutral round trip, and the named refusals.

## NC1 angle cutters (Claude, 2026-10-31)

Bob's aa5 (L4X4X5/16, NC1): GCmaker's EndCondition_NC1VM rules for Angle stock are reproduced by the GC Kernel
view's recipe (C#): an As Imported AK region crossing the root fillet's band (Through Fillet, automatic or set) is
drilled the fillet's radius past the leg, with the fixed 0.001 in overshoot; an edge-only preparation is its
welding wedge (cross-section between the waste overlap and the bevel, DSTV land from the back of the plate, or the
adjoining leg's source miter plane, swept along the edge). The wedge's bevel and the earlier triangle's face are the
same plane stored two ways, which the transverse chart now takes as one graph (above). tests/survey_gaps.rs checks
both cuts against the section's first moment below its top face (every point with y > L + z removed). C-Channel
NC1 cutters stay refused by name: GCmaker's channel rules (source-aligned end planes, straight-stock extension,
junction inflation, prepared flanges) are not reproduced, and the transverse chart refuses the first web cut of a
tapered-flange channel ("an end-cut surface patch is below the admitted partition resolution").

## Blind pockets in prism stock (Claude, 2026-10-31)

`pocket.rs` (inside the section domain): a right-prism tool along the body's prism axis that leaves through one cap
and ends inside the stock (a sketch cut part way into a Flat). Selected before construction unless both operands
are cardinal boxes (the rectilinear domain cuts those pockets). The result is exactly the stack A over [lo, d] and
A - B over [d, hi]: the bottom cap A, the top cap A - B, the floor A n B, one wall per segment of A stepping from hi
down to d where B covers it (no seam across a wall), one pocket wall per piece of B's boundary inside A. The 2D
arrangement is the section domain's (closed-form meetings, coincident pieces by direction, midpoints by winding);
loops chain by the least clockwise turn and a pinch is refused ("touches a boundary at a point"). A pocket from the
bottom is cut turned half a turn about the section's first axis (an exact rotation) and turned back. The volume is
checked against area(A) (d - lo) + area(A - B) (hi - d). Tests (tests/pocket.rs): rounded-rectangle, round,
triangular, rebate, flush, groove-across (top in two pieces), whole-top, counterbore about a through hole, bottom
pockets, through holes after a pocket (passage domain), all against area times depth to 1e-12; op 13 (the flat
nester's prism profile) still refuses a pocketed plate. Limits: one pocket per body (a pocketed body is no longer a
right prism; a second pocket is refused by name); op 12 history is not recorded for pockets (refused by name).

## Holes through conical walls (Claude, 2026-10-31)

`cone_bore.rs`: bodies with conical faces (tapered Pole stock, solid or hollow) take only this domain; its admission
names every refusal (conical tools stay refused). A full round cutter square to the stock's axis whose axial window
[z - rho, z + rho] lies strictly inside one outer wall face (cone or cylinder about the axis) and, hollow, one inner
wall face; starting outside the outer wall, ending in the bore (one wall), past the far side (through), or in solid
stock inside the material deeper than every entry point (blind: the bottom disk joins the hole's wall). On a cone
the loop is the new `Curve::ConePair` (wire type 8): per angle t round the hole, the height h along its axis solves
A h^2 + 2 B h + C = 0 with A = |a'|^2 - k^2 (a.b)^2, B = D'.a' - k R0 (a.b), C = |D'|^2 - R0^2 (D = Q(t) - p, primes
across the stock axis b, R0 = r + k D.b): not a conic. Its discriminant is a trigonometric polynomial of degree two
(exact from eight samples, proved positive all the way round). On a cylinder it is the CylinderPair. Clearances: a
loop's points lie within rho of the hole axis along b exactly (h a has no part along b when a is square to b) and
within asin((offset + rho) / R_min) of the hole's direction about the axis, so rings, seams and earlier holes are
cleared by those windows; other faces by boxes (pair edges: across both axes within rho, along the hole at
sqrt(R^2 - w^2) for the wall's radius range). Mass on the cone: Green's theorem in (theta, s) as for its parallels,
G(s) g(theta) dtheta along the curve by quadrature. Tests (tests/cone_bore.rs) check one wall, both walls off
centre and turned, solid through and blind, three holes in turn and one on the straight base against independent
disk quadratures of the walls' depth (1e-9). Not done: a hole on the revolve's seam (the seam is shared with the
adjoining section and not moved; refused by name), holes not square to the axis, holes ending inside a wall.

## Cut history in more domains; op 26 recorded cut (Claude, 2026-11-02)

Persistent names on every stock (Bob, contract in GCm GCmaker/Handoff/Claude-GcKernel-PersistentReferences.md):
op 12 history is now recorded, from the construction itself (`trace.rs`), in five more domains:
- **radial** (pipe holes) and **cone_bore** (holes through tapered or straight Pole walls, through cuts): the body is
  kept; each wall crossed gains one closed meeting loop [wall, cutter cylinder, n] (n counts that wall's crossings
  along the cutter) whose one vertex is the new role **42 MEETING_POINT [i, j, n] k=3**; each wall passage a seam
  TOOL_SEAM [entry wall, j] and a face from the cutter's cylinder. No history (the cut still runs) when the radial
  bore moves a wall's seam, or for a blind cone bore.
- **capsule** (slots): through a box, the box is rebuilt by op 10's extrusion with the slot as a hole, so its records
  are op 10's mapped back (sides, their cap edges, rulings and corners unchanged; caps modified; the hole from the
  tool); through anything else, each wall's cavity is spliced (`passage::splice_cavity`): the body kept, its two
  wall faces modified, cavity sides FROM_TOOL_FACE, cap copies MEETING, rulings FROM_TOOL_EDGE, their ends
  TOOL_EDGE_END. No history for the transverse cross-cut route.
- **box_cut** (rectangles and polygons): the passage splices as above; none for the transverse route or when coplanar
  cap patches had to be joined first (several parents of one face).
- **Selection change:** a plain cardinal box cutter making a clear passage through the planar walls of a body with
  no curved faces (`box_cut::clear_box_passage`, decided before construction) now runs in box_cut's passage branch
  instead of the rectilinear cell engine, whose cell seams split every face a cut touches. Same material removed
  (tests check exact volumes); every other rectilinear cut is unchanged. This changes op 3's body for those cuts
  (cleaner topology), the only op 3 change.
Op 12's refusals are unchanged in kind; a domain that records history now computes the body once and then reports
whether this cut's history was recorded.

**Op 26 CUT_RECORDED:** op 3's request; reply = op 3's body, then `1` and op 12's envelope (records and outcomes)
when the one selected domain recorded this cut's history, otherwise `0, n` and the n UTF-8 bytes of the reason. The
cut runs once; nothing is retried, matched or guessed. The host keeps names only while every cut carries records.
Tests (tests/provenance.rs): first interior rectangle, slot and turned square in a plate; hole, slot, rectangle and
square in sequence through an angle's leg; hole then rectangle through both walls of a tube; radial bores through one
and both pipe walls and a bore through both conical walls of a tapered tube; op 26's reply against op 3 and op 12.
Every case through the shared checker (each record against its geometry, generated names unique, outcomes re-derived)
with op 12's body byte-identical to op 3's.

## Union, ops 24 and 25 (Claude, 2026-11-02)

The public entries are `union(body, tool, tolerance)` (op 24: two bodies back to back as op 3 -> the union, one solid)
and `union_named` (op 25: op 24's reply, then op 12's provenance envelope with its outcome section; roles 30-41 with
body = the first operand and tool = the second). Same philosophy as the cut (`unite.rs`): one exact domain chosen from
the operands before anything is built, its failure returned, never retried; everything else refused by name.

### Boss on a planar face (`unite_face.rs`; history recorded)

The SolidWorks boss: the tool's faces lying on the plane of one body face F, facing it, form the base (shared seams
cancelled: a revolved disc's halves). Admitted when (1) the base is one region (outline plus holes) strictly inside F:
no closed-form meeting of any base edge with any edge of F within the tolerance, the outline inside F's region, no loop
of F inside the outline; (2) the boss lies above F's plane and meets it only in the base: used vertices at height
>= -tol, those within tol base vertices; other edges (lines, arcs, ellipse arcs) without interior height extremes at or
below tol; other faces planes off F's plane, or cylinders, elliptic cylinders, cones, spheres or tori about an axis
square to F (heights extreme on their edges); (3) every other body face's conservative analytic box in F's frame
(planes: their edges; cylinders/cones: axis segment widened by the radius; spheres/tori: whole) is at or below the
plane, above the boss's highest point or beside the boss's box. Any body: plates, webs, flanges, legs, sloping faces.
The body is kept entity for entity; F gains the base outline as an inner loop (the base's own edges); base holes become
new faces on F's surface; the tool's other faces/edges/vertices are added unchanged. Volume checked (body + tool).
History: every body entity FROM_BODY_* [i, 0] except F [F, 1] (and its islands, F split); tool faces other than the
base FROM_TOOL_FACE [j, 0]; every kept tool edge/vertex FROM_TOOL_*; base faces and seams consumed.
Refusals: "the boss's base touches or crosses the edge of the face it stands on, or a hole in it", "the boss's base is
not on the face it is sketched on", "a hole or edge of the face lies under the boss's base", "the boss reaches below
the face it stands on (a boss into the body is not this domain)", "the boss touches the face's plane outside its base",
"the boss dips to or below the face's plane along an edge", "the boss has a <kind> face whose height above the face it
stands on is not bounded by its edges", "another part of the body comes within the boss's room above the face (the boss
would run into it)", "the body has a <kind> face; its clearance from the boss cannot be proved yet", "the boss rests on
more than one plane of the body", "the boss's base is not one region", "the boss's base meets itself at a point". The
room test is conservative (boxes): it can refuse a boss that would fit.

### Layered right prisms (`unite_stack.rs`; history not recorded: op 25 refuses by name)

Both operands layered along one cardinal axis k: every face a cap (plane square to k) or a wall (plane, cylinder or
elliptic cylinder parallel to k) whose edges run square to k or along it. Each layer's section is read exactly from the
walls' vertical edges crossing its middle (paired along each wall's support, material left of k x outward), the stack
proved by layer areas times thicknesses equal to the measured volume. The union arranges every slab of both operands
(the section domain's split: closed-form meetings, points within tol one vertex, coincident pieces one span, other
slabs classified by midpoint winding), stands walls where material lies on one side, builds up/down level faces from
the material below and above, joins walls of one support and facing across shared edges (spikes removed, splits made
holes, one seam kept on bands round a periodic support) and joins edges of one support through vertices only they use
(`unite_tidy.rs`). Surfaces are the operands' own faces. Covers overhangs, flush bosses (one L-shaped side face), bars
through plates, overlapping/stacked bosses, bosses over holes and pockets, I-beam blocks along the beam. Volume checked
against the layer areas. Refusals: "the boss does not reach the body: their union would be two separate solids (one
solid only)", "the bodies would touch only along an edge (non-manifold result)", "the union would touch itself at a
point of its outline (non-manifold result)", "the union would leave separate solids or enclose a void (one solid
without voids only)", the section domain's ellipse-meeting refusal.

### Overhang on a non-layered body (`unite_overhang.rs`; history recorded)

Tried after the face and stack domains decline: a boss whose base (as the face domain's) rests on exactly one planar
face F and crosses F's outline transversally, on any body (a round boss over an I-beam flange tip or corner). Both
outlines are split in F's plane by the section arrangement; flush pieces (the base running along F's edge) are refused;
the boss must lie above the plane with nothing of the body in its room (the face domain's checks), and under the
overhang nothing of the body may reach the plane: no other face on the plane near the boss, and every other face
reaching the plane near it a plane or cylinder whose edges reach the plane only at F's own vertices (closed-form
extremes); no base corner on F's outline. Result: F less the base (a face per region), the overhang as the boss's
underside on the base's plane, crossed edges split at the crossings (pieces lifted from the arrangement, lines between
their vertices, arcs/ellipse arcs on their supports) and given to the faces beside them; everything else of both
bodies unchanged. Volume checked (body + tool). Refusals: "the boss's edge runs along the face's edge (flush) on a body
that is not a layered prism stack; not supported yet", "the boss rests on more than one face of the body", "a curved
part of the body reaches the plane under the boss's overhang", "another part of the body reaches the plane under the
boss's overhang", "a corner of the boss's base lies on the face's outline", plus the face domain's boss and room checks.
History (op 25): kept entities [i, 0]; F's regions [F, 1] (split when several); faces beside a split edge [i, 1]; the
underside FROM_TOOL_FACE [base face, 1] (op 25 refuses when the base is several faces); split edges' pieces [e, 1]; each
crossing vertex the new role 43 EDGE_CROSSING [body edge, tool edge, n], n counting crossings of that pair along the
tool edge's own direction (generated names are checked unique).

Not supported yet: unions none of the domains take (bosses on cylindrical faces such as a pipe's; a boss sunk partly
into a non-layered body; flush bosses on non-layered bodies), history for the stack domain, spline and
Bezier operands. Tests: `tests/union.rs` (closed-form volumes, topology check and tessellation of
every result, entity-for-entity preservation in the face domain, op 24/25 through the protocol, a hole cut after a
boss).

## Bores across cap seams (Claude 2026-11-03)
A transverse cut (endcut_cross) leaves the faces it crosses split into cells on one plane. A later bore whose corridor
meets such a seam on its caps joins the axis-normal faces it reaches, grouped by plane and outward side, across the
edges they share (`passage::join_cap_seams`, unite_tidy's exact topology edits; no point moved). Decided from the
operands before admission; the joined body is admitted as any other. Such a bore records no history yet (op 26 reason
"... across cap seams joined by this cut ..."). The transverse engine itself keeps its cells: joining them globally
broke later transverse cuts (patch coverage), so the join is local to the bore. Regression: tests/cap_seams.rs (Bob's
Angle Corner Cutout).

## Near-cardinal directions are cardinal (Claude, 2026-11-05)

Declared normalization contract, under Bob's rule of 2026-10-28 ("what agrees within the kernel tolerance is the
same, and the status says so"). Stock built by rotated, tapered or filleted profile math carries plane normals and
axes off the coordinate axes in the last bits: a C12x20.7 flange end face's (-1, 1.3e-16, 0), a C6x13 inner web's
(-2.6e-17, 0, -1). Two places now take a direction within 1e-12 per component of a coordinate axis as that axis.
Points never move; over any admitted part the change lies far inside the tolerance.

- **Transverse end cuts** (`endcut::snap_directions`, in `planar_cross_cut`, `cross::cut` and `cross::cut_along`):
  the operands' plane normals and x axes and cylinder axes are snapped inside the domain, after selection, so
  routing is unchanged. The captured pipe and rect-tube 45-degree miters, refused before, now cut. The rat holes'
  plane-fillet crossings come out as circles (arcs) instead of near-circular conics.
- **Round passages** (`passage::cap_normal`): such a face is a cap. A planar side face inside the corridor's box is
  excluded when it lies provably farther than the radius from the corridor's axis segment.

Fixtures: `parts/channel_end_partition` (Arc Steel B2025's first AK end cut, refused on b537031 as "an end-cut
surface patch is below the admitted partition resolution") and `parts/channel_web_hole` (Valley Industrial B24's
web hole).

## Stepped transverse faces (Claude, 2026-11-05)

The transverse coverage check (`verify_interfaces`) took every face parallel to the partition planes as a rectangle.
After one end cutter, a channel's end face can be an L or stepped outline where two material sections meet at one
station. Such a face is now any single loop of edges parallel to the y and z axes (`Outline`), with positive
oriented area. Its membership is decided even-odd on a grid that holds all its vertex coordinates. The coverage
equation (left sections minus right sections equals the signed patches, cell by cell) is unchanged. Fixture:
`parts/channel_stepped_end` (Bob Cody 12004SL0's second AK cutter).
