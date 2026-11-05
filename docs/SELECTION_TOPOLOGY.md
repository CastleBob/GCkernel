# Read-only selection topology

Owner: GC Kernel, selection module. Approved by Bob on 2026-11-05 to fix fragmented
face and edge measurement in the Rust viewport. This catalog does not modify the
B-rep, Boolean operations, save format, protocol, or tessellation.

`selection::build` accepts a body and positive finite validation tolerance. It
validates serialization and topology before indexing. The result contains face
regions, ordered edge chains with their analytic curves, and raw-ID-to-region
maps. IDs belong to one body revision; clients must retire them on replacement.

Face regions require shared-edge connectivity and identical oriented analytic
supports under the existing display seam predicate. Every member is compared
with the region representative to prevent accumulated drift. Disconnected
coplanar faces remain separate. Holes and distinct tangent supports remain.

Boundary chains require degree-two connectivity, the same two adjacent face
regions, and one line or circle support. Traversal checks endpoints, direction,
and circular sweep. The roundoff allowance is 64 machine epsilons, independent
of display chord or modelling tolerance. Internal same-support divisions are
not selectable edges. Other curve types remain individual edges; no spline
concatenation, approximate healing, or guessed recovery is provided.

Consumers highlight every raw member, measure the whole analytic chain, and sum
the member face areas (exact only when every member integral is exact). Point
targets use chain endpoints and whole-edge midpoints; periodic seam vertices
are not corners. Invalid catalogs return explicit errors and leave the body
unchanged. This is selection identity, not persistent topological naming.

Focused regressions cover unchanged boxes, split straight edges, fragmented
circles, nearby distinct supports, disconnected coplanar faces, invalid input,
and the rebuilt Angle part. The Angle fixture has 83 raw faces and 171 raw edges;
its catalog has 16 face regions and 41 boundary chains. Its top plane and slot
semicircles are checked analytically, including area and arc sweep.

Fixture provenance: tests/fixtures/selection/angle.f64 is our GCcore rebuild of
GCm GCmaker/TestParts/Angle.gcmprt captured 2026-11-05, little-endian neutral f64
values. SHA256: F2C14444A3671148DB7D644C9964CE27CFA6419BADA50B8955DDBA19507612E9.
No third-party implementation source or dependencies were used.
