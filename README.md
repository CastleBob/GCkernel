# GC Kernel

GCmaker's own solid-modelling kernel, in Rust (Bob, 2026-10-27). Clean room: our own code from public
mathematics; no third-party geometry source read, copied or translated. No dependencies.

Ownership: foundation (geom, frame, neutral, topology, extrude, mass, rebuild, tessellate, protocol, worker,
tests/fixtures) Claude; `src/boolean/` written by Astra, maintained by Claude since 2026-10-28 (Bob put Claude on
the modeler and cuts, Astra on the GC Kernel viewport). Plan and coordination: GCm `GCmaker/Handoff/Claude-GcKernel-BattlePlan.md`,
`Claude-GcKernel-Status.md`, `Claude-GcKernel-FilletCuts-Request.md`.

## Build and test

Independent STEP import (2026-10-31): `step::read_file` and worker op 16 read
supported standalone analytic solids directly into our Body, in inches. See
[scope, wire contract and verification](docs/STEP_IMPORT.md). General splines,
assemblies and IGES remain unfinished; this does not switch application intake.

```
cargo test --release
cargo build --release      # target/release/gckernel-worker.exe
```

GCmaker's developer-only lab (`GCmaker/Services/GcKernelLab.cs`) runs
`target/release/gckernel-worker.exe` unless `GCKERNEL_WORKER` names another. One request per process,
20 s deadline; nothing it returns enters a GCmaker document, saving or machining.

## Worker protocol v1 (little-endian, inches)

Request `GCKQ`, u32 version 1, u32 op, u32 units (1 = inch), f64 tolerance, u32 count, count f64.
Reply `GCKR`, u32 version 1, u32 status (0 ok, 1 unsupported, 2 invalid input, 3 internal), u32 count,
then count f64 (ok) or count UTF-8 message bytes.

Ops: 1 check, 2 extrude, 3 cut, 4 measure ([volume, area, 1 closed form | 0 quadrature]),
5 rebuild prism, 6 tessellate ([chord, body] -> display mesh and edge polylines, `src/tessellate.rs`; display only),
7 transform (one rigid motion and a batch of items -> the items moved; `src/frame.rs`, `protocol::transform`):
`[1, mode (0 T, 1 T inverse), T (m row-major 9, t 3), N, N x [kind, len, data]]` -> `[1, N, N x [kind, len, data]]`;
kinds 0 point, 1 vector, 2 frame (o x y), 3 rigid (reply: T after it), 4 curve, 5 body, 6 frame-local point
(frame, x, y; reply kind 0), 7 surface. One bad item fails the batch with its index.
Curves: 0 line, 1 arc, 2 general B-spline (refused), 3 conic, 4 cylinder pair,
5 polynomial cubic Bezier (four xyz points), 6 rational cubic Bezier (four xyz
points then four positive weights); the Bezier parameter domain is [0, 1].
Surfaces: 0 plane, 1 cylinder. See `src/protocol.rs`, `src/neutral.rs`, `src/geom.rs`.

Display polylines omit internal same-plane/same-cylinder divisions and periodic seams. Their IDs remain the
original body edge IDs (possibly sparse); exact topology and mesh boundary samples are retained. Distinct
supports, including tangent fillet transitions, remain visible. Multi-hole face triangulation chooses bridge
endpoints using the local interior sector when earlier holes repeat a vertex. These are display-only rules.

## Precision

Floating point with explicit tolerances and whole-interval mathematical bounds where stated; not exact
arithmetic and not interval-certified. Valid input outside a supported domain is refused (Unsupported),
never approximated, repaired or retried another way. `mass::measure` integrates whatever body it is
given; every protocol path checks the body first.

Geometry primitives (`src/frame.rs`): Vec2 / Point2, Point3 beside Vec3 (explicit conversions; bodies keep Vec3),
Frame (right-handed, z = x cross y) and Rigid (rotation + translation; points move, vectors only turn). Supplied
unit lengths, right angles and rotation matrices are admitted within `frame::ADMIT` (1e-12) and otherwise refused,
never normalized; no mirror, scale or shear is representable. Angles in radians.

## Precision rule (Bob 2026-10-28)

What agrees within the request tolerance is the same, and the status says so (for example a hole whose wall
touches a fillet exactly at its tangent line: `geom::pair_trim_touches`). Dimensions are never nudged to dodge a case.

## Fixtures

`tests/fixtures/` holds bodies captured from GCmaker exactly as sent (neutral f64 arrays + manifest):
the I-beam rat-hole part, circle over fillet, a shop-cut survey (tube and pipe), and Bob's real parts
(STR-1, Rect Tube) with per-step Eyeshot references. Regenerate from GCm `GcKernelLabTests` with the
`GCKERNEL_FIXTURES`, `GCKERNEL_SURVEY`, `GCKERNEL_PARTS` and `GCKERNEL_PARTS_REF` variables.
