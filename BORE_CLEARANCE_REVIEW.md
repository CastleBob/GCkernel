# Codex review of b537031 round-passage clearance

This branch contains two intentionally failing regression tests only. No production
algorithm change; do not merge/promote it as a verified fix.

At a 1e-6 tolerance, the axis segment (0,0,0)-(0,100000,0) and edge
(-0.0001,0,0)-(0.0001,100000,0) cross at (0,50000,0). segment_distance
returns 0.0001 because a*e-b*b cancels to zero and takes the parallel branch.

A planar CCW strip with vertices (-0.0001,0,0),(-0.00009,0,0),
(0.00011,100000,0),(0.0001,100000,0) contains the axis point (0,47500,0).
The test confirms inside() is true there. Nevertheless clear_of_corridor returns
true for the Y-axis segment and radius 1e-5 at tolerance 1e-6. Thus it is not a
conservative exclusion. This is a helper-level admission counterexample, not a
claim of a known corrupted customer part.

Run: cargo test --release --lib codex_clearance_review -- --nocapture
Observed: both tests fail on b537031. Log in Codex TEMP:
GCm-b537031-clearance-review.log. Initial fixture compile used nonexistent Vec3
axis constants; corrected to explicit vectors before the reproduced failures.

Claude owns the active boolean lane. Please correct the exclusion to remain
conservative under cancellation and retain this crossed-face refusal alongside
the newly admitted channel_web_hole fixture. GCm's installed workers stay on
f21974dc until the new pair passes review. No source tolerance change, geometric
repair, external algorithm implementation or dependency was used in this review.
