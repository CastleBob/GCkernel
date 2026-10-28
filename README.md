# GC Kernel

GCmaker's own solid-modelling kernel, in Rust (Bob, 2026-10-27). Clean room: our own code from public
mathematics; no third-party geometry source read, copied or translated. No dependencies.

Ownership: foundation (geom, neutral, topology, extrude, mass, rebuild, tessellate, protocol, worker,
tests/fixtures) Claude; `src/boolean/` written by Astra, maintained by Claude since 2026-10-28 (Bob put Claude on
the modeler and cuts, Astra on the GC Kernel viewport). Plan and coordination: GCm `GCmaker/Handoff/Claude-GcKernel-BattlePlan.md`,
`Claude-GcKernel-Status.md`, `Claude-GcKernel-FilletCuts-Request.md`.

## Build and test

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
5 rebuild prism, 6 tessellate ([chord, body] -> display mesh and edge polylines, `src/tessellate.rs`; display only). Curves: 0 line, 1 arc, 2 B-spline (refused), 3 conic, 4 cylinder pair.
Surfaces: 0 plane, 1 cylinder. See `src/protocol.rs`, `src/neutral.rs`, `src/geom.rs`.

## Precision

Floating point with explicit tolerances and whole-interval mathematical bounds where stated; not exact
arithmetic and not interval-certified. Valid input outside a supported domain is refused (Unsupported),
never approximated, repaired or retried another way. `mass::measure` integrates whatever body it is
given; every protocol path checks the body first.

## Precision rule (Bob 2026-10-28)

What agrees within the request tolerance is the same, and the status says so (for example a hole whose wall
touches a fillet exactly at its tangent line: `geom::pair_trim_touches`). Dimensions are never nudged to dodge a case.

## Fixtures

`tests/fixtures/` holds bodies captured from GCmaker exactly as sent (neutral f64 arrays + manifest):
the I-beam rat-hole part, circle over fillet, a shop-cut survey (tube and pipe), and Bob's real parts
(STR-1, Rect Tube) with per-step Eyeshot references. Regenerate from GCm `GcKernelLabTests` with the
`GCKERNEL_FIXTURES`, `GCKERNEL_SURVEY`, `GCKERNEL_PARTS` and `GCKERNEL_PARTS_REF` variables.
