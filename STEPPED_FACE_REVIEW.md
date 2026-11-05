
## Codex: stepped transverse face review scope

2026-11-04: Reviewing kernel f939898 on an isolated codex/stepped-face-review
branch, test code only in endcut_cross.rs. Exact rectilinear L-prism controls
will check the new coverage admission against an analytic material section,
including a filled-in notch, duplicate/missing cap and reversed orientation.
No production source edits or active checkout changes. The normalization review
remains unresolved; exact-axis tests do not clear that separate issue.

## Result

Two independent analytic tests pass on f939898 and both fail against the previous
3755a53 rectangular-only coverage function. The original code was restored after
that comparison; no production implementation changes are included here.
A 3x3 L section missing its 2x2 corner has area 5; length 2 gives volume 10.
Trimming length 0.5 leaves volume 7.5 and closed genus-zero topology. Tests run
at scales 0.01, 1 and 25.4, with scaled request tolerance. Direct coverage tests
reject the filled-in notch, duplicate cap, missing cap and reversed orientation.
The cut test enters cut_local explicitly, independently of routing/normalization.

Four existing channel/stepped-channel/miter checks also pass. Logs:
TEMP/GCm-stepped-face-before.log (2 red), GCm-stepped-face-review.log (2 green),
GCm-stepped-face-regressions.log (4 green).

Review conclusion: the exact rectilinear stepped-face extension has focused
independent evidence and preserves these coverage refusals. The separately
reported 3755a53 direction-normalization failures remain unresolved; this does
not approve the entire current kernel or prove all tolerance-degenerate outlines.
No main checkout, shared worker, installed pin or source geometry was changed.
Claude may take this test-only commit after reviewing it.
