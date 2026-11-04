# Codex review of b537031 round-passage clearance

Review branch: codex/bore-clearance-review. The original commit 4255a27 records
both failures below. The follow-up adds a conservative refusal guard; this is a
review delivery, not a main merge or an installed-worker promotion.

At tolerance 1e-6, axis segment (0,0,0)-(0,100000,0) and edge
(-0.0001,0,0)-(0.0001,100000,0) intersect at (0,50000,0). On b537031,
segment_distance returns 0.0001 because a*e-b*b cancels to zero.
A planar CCW strip with vertices (-0.0001,0,0),(-0.00009,0,0),
(0.00011,100000,0),(0.0001,100000,0) contains axis point (0,47500,0),
confirmed independently by inside(). clear_of_corridor nevertheless excludes
it for radius 1e-5. This is a helper-level counterexample, not a claim about a
known corrupted customer part.

The correction makes the private distance helper return Option<f64>. A
nonparallel determinant within 16 machine epsilons of its product scale supplies
no clearance evidence. Exact cardinal parallelism is recognized from zero
components, not from a rounded determinant. Nonfinite arithmetic also supplies
no evidence. The corridor caller refuses exclusion on None, preserving its
existing unsupported-geometry outcome. No guessed distance, alternate algorithm,
geometry repair, tolerance increase or external dependency is introduced.

Verification (release, locked):
- cargo test --lib codex_clearance_review: 3 passed, including crossed-face
  refusal, cancellation refusal, exact parallel and well-conditioned distances.
- cargo test --test bore_clearance --test cap_seams --test pocket --test provenance:
  29 passed (1 + 2 + 7 + 19). The channel web-hole acceptance remains passing.
- Original 4255a27: both counterexamples fail, log
  TEMP/GCm-b537031-clearance-review.log.
- Follow-up logs: TEMP/GCm-bore-clearance-guard-{unit,regressions}.log.

Research: https://www.cs.cmu.edu/~quake/robust.html explains why near-zero
floating-point determinants require explicit uncertainty handling. Only the
public mathematical overview was consulted; no linked implementation was read
or reused. This local guard addresses the demonstrated cancellation class; it
is not a proof that every distance/clearance calculation is numerically certified.

Claude owns the active boolean lane. Please review/cherry-pick both commits
(or their combined change) and retain the crossed-face refusal alongside the
channel acceptance. The active kernel checkout was untouched. GCm's installed
pair remains GCcore 9f009e9 / GCKernel f21974dc pending reviewed integration.
