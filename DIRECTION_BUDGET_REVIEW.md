# Direction normalization: existing cylinder residual budget

Review of a08814a, test-only codex/direction-budget-review.
The three prior plane counterexamples now pass. Two new cylinder cases fail:
- length 90 in, tol 1e-10: returned edge/support departure 1.3500001116995008e-10;
- length 900000 in, tol 1e-6: departure 1.3500000000110646e-6.

Both originals pass endcut::validate and topology::check. The circle boundary
has radius 1 + 0.9*tol while its cylinder support has radius 1, an already
admitted residual. A Y rotation of 5e-13 leaves the original valid. The new
axial-reach bound permits another 0.45*tol when snapping that axis; combined
residual is 1.35*tol. Half the tolerance is not automatically unused budget.

Required correction: certify the whole boundary against the proposed cylinder,
or bound both its existing support residual and normalization departure within
the total requested tolerance. Do not increase tolerance or repair the source.
The radius contribution to axis-rotation bounds also warrants a mathematical
check; the current reach-only product is not a general whole-support proof.
This review has no production fix and does not claim an invalid complete Boolean
was accepted or that a customer file was corrupted. These are helper regressions.

Run: cargo test --release --lib codex_direction_snap_review
Result: 3 existing tests pass, 2 new tests fail. Full log:
TEMP/GCm-direction-budget-review.log. Active checkout and installed workers untouched.
## Isolated correction candidate - 2026-11-04

Branch codex/direction-budget-fix keeps the original red cases above intact.
Cylinder support motion now includes radius: for unit axes a/b, the difference
of axial projectors aa^T and bb^T has norm at most 2*|a-b|. Every point on the
bounded cylindrical patch is at most axial_reach+radius from its support origin,
so the resulting radial-distance change is bounded by their product. The
existing numerical margin is retained; this bound must fit half the request
length tolerance.

A separate certificate checks every whole edge curve against the PROPOSED
cylinder with the kernel's existing analytic topology::off_surface predicate,
and checks every incident vertex. Both must fit half tolerance. This includes
the input's existing residual; no unused budget is assumed. Points, curves,
indices and loops remain unchanged. Uncertified candidates retain their original
support under the existing normalization policy. The emitted note says
"largest certified departure bound", and uses the accepted half-tolerance bound
for cylinders rather than misreporting the additional axis motion as total error.

Seven focused checks pass: original three planes, original two cylinder failures,
a length-40 variant that specifically needs the new boundary check even though
the motion bound fits, and a positive radius-residual cylinder that still snaps,
reports and preserves its boundary. Five affected existing checks pass: channel
end partition, stepped channel, bore clearance and both miter_pair cases.
Logs: TEMP/GCm-direction-budget-fix.log and -regressions.log.

No new external source/dependency or numerical algorithm was imported. The
projection norm uses the triangle inequality; whole-curve predicates are the
kernel's existing independently authored topology checks. This is a candidate
correction, not Claude acknowledgment, a complete Boolean proof or a product
promotion. Active GCKernel and installed workers are untouched. Worker/GCm
interop against this corrected dependency is the next separate check.
