# Codex review of b537031 round-passage clearance

Branch: codex/bore-clearance-review. Review delivery only; active master and the
installed GCm pair remain unchanged. Claude owns the boolean lane.

## Reproduced failures and superseded guard

Commit 4255a27 preserves two failures on b537031. Axis (0,0,0)-(0,100000,0)
and edge (-0.0001,0,0)-(0.0001,100000,0) intersect, yet the normal-equation
distance is 0.0001 because a*e-b*b rounds to zero. An independently checked
planar strip crossed by the axis is then wrongly excluded from the bore.

Commit c5ef911's denominator guard fixes those two cases but is INSUFFICIENT.
Commit 3ed94d9 preserves two further failing tests: with x endpoints +-0.01,
the nonzero determinant yields distance 5.763319672131423e-6 for an intersection.
A crossed strip is wrongly excluded at radius/tolerance 1e-6 when its long edges
are stored in +Y (the last loop occurrence is reversed). Do not adopt c5ef911 or
the red-only 3ed94d9 as the final correction.

## Current correction: verify the proposed separation

segment_clearance_lower_bound returns a lower bound, not the distance estimate.
The existing normal equations propose direction n only. For each original
endpoint, every dot-product multiplication and addition is rounded outward with
next_down/next_up. Endpoint intervals therefore enclose every projection along
each segment. The positive gap between those intervals is rounded downward.
An upper bound on |n| is computed by rounding squares, sums and square root up;
the final gap/norm quotient is rounded down.

Proof: if every q projection exceeds every p projection by g, then
n dot (q-p) >= g. Cauchy-Schwarz gives |q-p| >= g/|n|. Endpoint projection
extrema bound the entire line segment, so this applies to all point pairs.
Using a lower gap and upper norm preserves that lower bound. The symmetric
ordering is checked too. An inaccurate proposed direction can weaken the proof
but cannot establish separation of intersecting segments. Overlap returns zero;
uncertainty or nonfinite bounds refuses exclusion. Existing near-parallel refusal
is retained. No guessed distance, alternate geometry/recovery path, tolerance
increase, source repair or dependency is introduced.

This certifies the segment-boundary evidence against represented f64 endpoints;
it is not a blanket certification of the rest of clear_of_corridor's plane/interior
logic or every boolean operation. These are helper-level counterexamples, not a
claim that a particular customer part was corrupted.

## Verification

- cargo test --locked --release --lib codex_clearance_review: 6 passed, including
  both pairs of counterexamples, positive parallel/skew bounds, and 72 analytic
  midpoint intersections across scales, axes and endpoint orderings.
- cargo test --locked --release --test bore_clearance --test cap_seams
  --test pocket --test provenance: 29 passed (1 + 2 + 7 + 19), including the new
  channel web-hole acceptance and existing cut history.
- Logs: TEMP/GCm-clearance-certificate-{unit,regressions}.log.
- Earlier red logs: GCm-b537031-clearance-review.log and
  GCm-clearance-nonzero-red.log. Original failures remain in git history.

References consulted as public math/API documentation only:
- https://www.cs.cmu.edu/~quake/robust.html (determinant roundoff)
- https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-curve-prop.html
  (convex-hull containment; the segment case is linear interpolation)
- https://doc.rust-lang.org/std/primitive.f64.html#method.next_up
- https://doc.rust-lang.org/std/primitive.f64.html#method.sqrt
  (adjacent representable floats and correctly rounded square root)
No linked implementation source was read or reused. Projection-bound reasoning
and code were independently written. Please review the complete branch change
before any main merge/promotion. GCm remains GCcore 9f009e9 / GCKernel f21974dc.
