# Annular pipe-window review candidate

Base: 04bbd40. Branch: codex/annular-window. Owner: Codex isolated candidate;
Claude retains active kernel ownership. No installed worker or GCm pin change.


## Annular pipe windows (Codex review candidate, 2026-11-04)
The AnnularWindow domain is selected for a two-cylinder coaxial sector tool
with a partial circular edge and a stock cylinder parallel to that axis.
Its admission proves a plain four-face pipe and a six-face annular sector
from every face, loop, edge and corner. The tool must overreach both pipe
wall radii and be strictly clear of both end caps. Axis lines, angular recipes
and cap planes agree exactly; there is no near-axis or near-plane substitution.
The angular interval and its complement must exceed the stated clearance
budget. Generic topology validation is followed by this stricter domain proof.

The result retains the analytic pipe material outside the sector: two walls
with window loops, two annular end caps and four planar window walls. A fresh
periodic parameter seam lies halfway around the angular complement. This is
the admitted plain-pipe boundary, not a regenerated stock with earlier features
discarded: any earlier opening, trim or extra boundary refuses before building.
Multiple windows, normalized/mitered ends, conical Pole stock, blind wall
pockets, full-turn severing and tangent/end-touching cuts are outside this domain.

The result has no ancestry records because the seam representation changes.
cut_recorded returns its actual missing-history reason; cut_named refuses before
construction. There is no retry through an unnamed operation. Other domains
and tolerances are unchanged. Inputs are immutable. Arithmetic uses f64, not
exact/interval-certified arithmetic; curves and surfaces stay analytic.

The independent volume is pi*(Ro^2-Ri^2)*L minus
theta*(Ro^2-Ri^2)*(z1-z0)/2. Surface area removes the two wall rectangles and
adds two annular sectors and two radial rectangles. These derive from public
cylindrical coordinates (MathWorld equation 17, r dr dtheta dz), written here
independently; no third-party kernel source/dependency.

Verification: five new integration tests, 94 existing Boolean unit tests,
11 normalized-end tests, 19 provenance tests and 11 survey-gap tests pass
(**140 checks**). New cases cover seven sweeps (0.01, 30, 90, 180, 270, 350,
359.99 degrees), both directions, independent volume and area, actual clock
angle, neutral round trip, immutable inputs and closed outward display meshes.
Translated cardinal frames and 0.1/1/25.4 scales pass. Both revolved and extruded
plain-pipe representations pass. Insufficient radial reach, touching/crossed
end caps, prior windows, inverted walls and a slightly displaced cap (inside
generic topology tolerance but outside the exact recipe) explicitly refuse.
Recorded/plain bodies agree; unavailable ancestry and strict-name refusal pass.

Logs: TEMP/GCKernel-annular-window-{regression,boolean}.log. No full-suite,
native EyeShot, saved SquareCutoutRevolved adapter, arbitrary-axis, persistent-name
or production-readiness claim. GCcore c5861a9 remains the original red reproducer;
a separately pinned consumer will verify that same test against this candidate.
