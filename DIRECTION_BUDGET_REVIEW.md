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