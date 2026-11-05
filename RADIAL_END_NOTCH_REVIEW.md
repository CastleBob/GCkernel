# Centered radial flat-end notch — Codex, 2026-11-05

Review branch codex/radial-end-notch from e0a0105. Primary ownership remains
Claude's Boolean subsystem; active checkout and installed pins are unchanged.
New radial_notch.rs and six-line domain dispatch addition implement the narrow
contract appended to src/boolean/README.md. This is independently written from
analytic cylinder intersection and existing project geometry representations.
No third-party implementation or dependency.

The near-end half disk lies at axial coordinate >=0; its intersection with each
pipe cylinder is one CylinderPair edge over t=0..pi. Split end-ring complements
and two radial-end connecting lines close the end cap and cutter cylinder.
Each new boundary is stitched by vertex identity, never proximity or welding.
A whole simple-shell proof excludes hidden prior openings and preserves distant
planes/cones without regenerating them. Exact frame/center/coaxial relations and
strict finite reach/clearance are required. No guessed recovery or tolerance edit.

Verification: six new focused tests pass, covering plain/distant-bevel stock,
both ends, three cardinal frames with translation, two physical scale factors,
independent half-disk volume quadrature, neutral round trips, outward closed
meshes, unchanged far faces/edges and unchanged operands. Seven refusal cases
cover previous/identical openings, bevel crossing, seam crossing, wall pocket,
excessive radius and offset center; inverted cutters are also rejected.
Op 26 computes the body with an explicit unavailable-history reason; strict
named op 12 refuses. Seven existing annular-window, five conical-bore and four
radial-bore regressions pass (22 distinct focused tests including the new six).
Logs TEMP/GCKernel-radial-notch-{final,regression,radial-regression}.log.

Original GCcore/GCm consumer checks are next, not yet claimed. In particular,
Tube_3in's repeated identical cutters are not admitted by this first domain;
Pipe Coped Chamfered's bevel-intersecting cutter also remains outside it.
Do not mark either complete. No installed promotion, native-parity waiver,
product publication or Claude acknowledgment.
