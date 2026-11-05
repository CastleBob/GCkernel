# Tangent straight/circular outer-rim chamfers

Codex candidate, 2026-11-04 night, branch codex/tangent-rim-chamfer from 174d77f.
This addresses the kernel blocker isolated by GCcore a3dbc1e. It does not add a
saved End Prep adapter, choose per-side ancestry or install a worker.

## Domain and construction

At a three-face vertex, exactly two selected convex chamfers share one planar
end face. One edge is straight between two planes; the other is circular between
the end plane and a cylinder. The remaining edge is an axial ruling. The skin
plane must be tangent to the cylinder with the same material side; end/skin
setbacks must meet at the same two proposed rail points. Only the outer rim
(radial setback toward the cylinder axis) is admitted. Fillets, unequal setbacks,
inner rims, non-tangent junctions and other mixed vertices refuse explicitly.

With cylinder section center C, radial direction rho and axial normal a, the
joint endpoints are T=C+(R+alpha*sp)*rho and Q=V+beta*sc*a. Holding angle fixed
in a cone's parametric equation gives a straight generator T-Q. The planar
chamfer's two rail endpoints are independently checked against T and Q before
they share those topology vertices. Residuals use at most tolerance/8; angular
admission is at most 1e-12 and tolerance/(8*incident reach). Existing complete
face, clearance, trim and topology checks remain mandatory. These are numerical
admission bounds, not interval arithmetic or an unlimited formal FP proof.

The original interference envelope for every circular rim was the entire
annulus. That conservatively rejected an inner boundary outside a rounded
corner's actual angular interval. For arcs at most pi radians, the same radial
and axial envelope is now intersected with two inward radial half-spaces. Their
normals are unit length; admitted arc-axis/skin-axis residual over the axial
range and tolerance expand the envelope conservatively. The existing subdivision
proof, budget and every tested input edge remain. Longer arcs keep the original
containing annulus. This changes only a containing clearance bound, never the
requested geometry or the failure algorithm.

Public mathematics reviewed: Eric Weisstein, [Cone](https://mathworld.wolfram.com/Cone.html),
parametric/implicit equations 1-6. The fixed-angle generator and sector half-space
derivations and implementation are ours; no third-party kernel implementation or
new dependency was used. A Rice lecture search result returned 404 and was not
relied on.

## Verification

All 39 tests in tests/blend.rs pass (35 existing, four added). The added capsule
stock cases cover both ends, unequal end/skin distances, reversed request order,
sequential opposite-end preparation, analytic volumes, neutral round trip and
closed outward display meshes. Unequal adjoining setbacks, mixed fillets,
non-tangent joins and a hole inside the chamfer sector remain refusals; input
bodies remain unchanged. Test log: TEMP/GCKernel-tangent-rim-regressions.log.

The unmodified GCcore a3dbc1e capability tests now both pass with this kernel:
24 pipe blend calls with volume/area/root-land checks and four rounded-rectangle
side/unit cases with independent offset-section integrals. The original test-only
branch remains red against 174d77f, preserving the reproducer. Consumer log:
TEMP/GCcore-tangent-rim-candidate.log. These checks do not claim arbitrary blend
support or complete native End Prep parity. Review/adoption is requested from
Claude; no acknowledgment or installed promotion is claimed.
