
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

## Exact input-plane arithmetic follow-up — 2026-11-05

Temporary logging at `endcut::hit` captured the two source planes in the
transverse frame. Logging was removed; runtime source has no diff. The
standard-library-only `tests/channel_support_arithmetic.py` preserves the
captured binary64 inputs and checks their mathematics with exact Fractions.
Three audit checks pass; the original Rust refusal check also passes.

This strengthens and qualifies the earlier construction-bound suggestion:

| Quantity | First plane | Second plane |
| --- | --- | --- |
| Exact-input intercept (decimal approximation) | -2.6146169744957527e-16 | -9.209377708898327e-33 |
| Stored computed intercept | -2.7755575615628914e-16 | 0 |
| Local intercept conversion error | -1.6094058706713886e-17 | 9.209377708898327e-33 |

The intercept conversion errors total less than 7% of the exact-input
intercept separation. The difference of the exact supplied supporting lines
is -2.6146169744957444e-16 at the band's low end and +2.1275767045623574e-17
at its high end. Their exact-input intersection is approximately
0.30664916816342186, inside the band. Thus these binary inputs do not define
one exactly coincident line; a bound on ONLY the local plane-to-line arithmetic
cannot establish that they were intended to be the same support.

This does not establish meaningful physical separation or prove the original
authored planes differ: earlier construction/frame transformations may have
introduced the discrepancies. Further work should trace those operations and
their error/provenance, or explicitly define and review a support-equivalence
contract. It must preserve genuinely distinct close surfaces and not infer
intent from small numerical separation alone. No tolerance floor, snapping,
fallback or geometry acceptance was added. This remains an unresolved cut.

Capture log: TEMP/GCKernel-channel-plane-audit.log. Script checks reproduce the
two logged graph values, prove the rational crossing, and quantify conversion
error. They audit the recorded values, not all current runtime constructions;
rerun the capture if that construction changes.

