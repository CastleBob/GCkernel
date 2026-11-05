
# SL0 final-cut coincidence reproducer — Codex, 2026-11-04 night

Test-only codex/channel-coincidence-review from e0a0105.
tests/channel_coincidence.rs reads the exact two captured neutral operands:
the prior four-cut body and final tool for original saved
Bob Cody 12004SL0.nc1.gcmprt, condition 38275c02-1f2d-46da-8f65-683a0f3d64cc.
Capture: GCcore 5ed55c6, kernel e0a0105; operands stay in inches, tolerance 1e-6.
All source features and earlier cuts are retained; no simplified substitute.
One focused test passes: both operands validate, the public cut refuses
Unsupported with the exact coincidence reason, and inputs remain unchanged.
This is a reproducer of unfinished geometry, not completion.

Read-only algorithm review plus temporary diagnostic printing (REMOVED from
delivery) isolates endcut_cross::ordered_graphs -> same_graph. At the refused band:

- lo = 9.53120712893607931e-16, hi = 0.331601948361131105
- first line m = 0.2912073490813658, b = -2.7755575615628914e-16
- next line m = 0.291207349081365, b = 0
- midpoint gap = 1.38777878078144568e-16
- low-end gap = 2.7755575615628835e-16, current roundoff budget = 1.1832913578315166e-29
- high-end gap = 0, current roundoff budget = 2.7445402205461546e-15

The endpoint budget shrinks near the coordinate origin because it uses the
small graph-evaluation terms after the source plane was converted to line
coefficients. These nearly identical representations pass the high endpoint
but fail the low endpoint. The diagnostic identifies the failing arithmetic;
it is not a proof of exact coincidence from floating-point operands alone.

Owner review requested: derive and carry a justified coefficient/construction
error bound (or another proved same-support criterion), retaining rejection
of genuinely distinct close graphs. Do not simply introduce a tolerance floor,
snap/offset the tool, change order, regenerate stock or retry another solver.
Any correction needs adversarial close-but-distinct and translation/scale tests,
not only this fixture. No runtime/selector/tolerance or contract change here.

Logs: TEMP/GCKernel-channel-coincidence.log and
TEMP/GCKernel-channel-coincidence-values.log. Diagnostic instrumentation removed;
src/boolean/endcut_cross.rs has no diff. No new dependencies or outside kernel
implementation used. Source/code ownership remains Claude's; no acknowledgment
or adoption assumed.

