# Review of near-cardinal normalization in 3755a53

Branch codex/direction-snap-review contains three intentionally failing regression
tests and this note only. No production implementation change; do not merge it as
a correction. Bob's existing within-tolerance equivalence policy is accepted.
The tests show that a fixed direction-component threshold is insufficient to
establish the required spatial bound.

Each test constructs an exact line-profile prism, rotates it by 5e-13 radians
about Z, and passes the existing endcut::validate at the requested tolerance
BEFORE snap_directions. The helper returns a body; topology::check then fails:

1. Length 1000 inches, tolerance 1e-10: boundary/surface departure 5e-10.
2. Length 1e8 inches, ordinary tolerance 1e-6: departure about 5e-5.
3. A ONE-INCH prism, tolerance 1e-6: move each plane origin 1e8 along its own
   tangent x direction (the same infinite plane, a valid parameterization).
   The input again validates. Snapping produces departure about 5e-5.

Thus even bounding the part's box alone would not cover a valid support origin
far from its patch. Points staying fixed does not mean their supporting surfaces
stay within tolerance. These tests disprove the helper's stated equivalence
claim; they do not claim that a complete public boolean accepted a corrupt output
or that a customer file was corrupted.

Run: cargo test --locked --release --lib codex_direction_snap_review -- --nocapture
Observed: all three fail after the original input validates. Log in Codex TEMP:
GCm-direction-snap-review-red.log.

Requested correction under Claude's boolean ownership:
- Use the request tolerance and certify spatial displacement over the bounded
  support/face domain before accepting a direction normalization. Explicitly
  refuse equivalences that cannot be established. Do not increase tolerance,
  repair input geometry or silently retry another path.
- Keep Bob's required reporting of which equivalences were used. Current
  snap_directions returns only a Body and no equivalence diagnostics.
- Align the documented domain with the actual implementation: snapped() clears
  small components of arbitrary directions, including oblique directions near
  coordinate planes, not only vectors near a coordinate axis.
- Retain the newly accepted channel and miter fixtures with independent geometry
  expectations, along with request-tolerance and support-origin regressions.

GCm's installed pair remains GCcore 9f009e9 / GCKernel f21974dc. Separately,
GCcore codex/current-worker-admission at 5fe85f5 pins projection-certified
GCkernel b0637e6; its known zero-GUID host admission blocker remains open.
This normalization review is not included in that combined review pair.
