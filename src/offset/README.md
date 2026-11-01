# Offset and shell (Claude, 2026-11-01)

Worker ops 20 SHELL and 21 OFFSET (layouts in `src/protocol.rs`). Clean room: built from the public mathematics below,
no third-party kernel source consulted or translated.

## What is exact, and how

- **Straight prisms of lines and arcs** (`prism.rs`): plates, plates with holes, rectangular tube, angle, channel and
  beam profiles. Inward offsets (rounded or sharp) and sharp outward offsets are prisms of the planar offset profile.
  The rounded outward offset is the true parallel body: side walls are the planar offset (convex profile corners
  become cylinders about the old rulings), each cap's band is the tube of radius d about the cap outline (quarter
  cylinders about lines, quarter tori about arcs, sphere pieces at convex corners), and at a concave line-line corner
  the two quarter cylinders meet along an exact ellipse in the bisecting plane. Volume in closed form:
  `L (A + P d + K d^2) + 2 (A d + P pi d^2 / 4 + 2 K d^3 / 3)` per closed cap band.
- **Full solids of revolution of lines and arcs** (`revolved.rs`): the planar offset of the whole section through the
  axis (section plus mirror image), cut back to r >= 0 and revolved exactly (`revolve.rs`). Volume by Pappus.
- **Shells**: closed (the inward or outward offset as a void, an inner shell of the same body: the neutral format and
  `topology::check` accept several shells), or open through removed faces: prism caps (a box opens through any face,
  being a prism three ways), or any revolved faces (both halves of a split surface), the cavity's wall ending flush
  on the removed face.

## Refused by name

Elliptic cylinders and general meridian surfaces (their offsets are not ellipses or in the kernel's surface set), arcs
whose radius the offset reduces to zero or below, sides eaten by their neighbours, loops that touch, cross or change
nesting after the offset (holes merging, a bore closing, a wall consumed), inward offsets or walls consuming the
thickness, concave sharp corners beside an arc in a rounded outward prism offset (the tubes meet in a non-conic
curve), removed faces that are not prism caps, half of a split revolved surface, removed faces meeting kept faces at
concave or tangent edges, part-turn revolves, and bodies outside both domains.

## Sources (public mathematics)

- Parallel bodies and Steiner's formula for the volume of the outer parallel body of a convex body (J. Steiner 1840;
  any convex-geometry text, e.g. R. Schneider, *Convex Bodies: The Brunn-Minkowski Theory*, ch. 4); the distance to a
  product set A x I is sqrt(dist_A^2 + dist_I^2).
- Offset (parallel) curves of lines and circles, and offset surfaces of planes, cylinders, cones, spheres and tori as
  surfaces of the same kind (standard differential geometry, e.g. do Carmo, *Differential Geometry of Curves and
  Surfaces*; Farouki, "The approximation of non-degenerate offset surfaces", CAGD 1986, for the exact cases).
- Rolling-ball (constant-radius) blends: pipe surfaces about lines and circles are cylinders and tori; two equal
  cylinders with intersecting axes meet in plane ellipses (elementary analytic geometry).
- Pappus's centroid theorem for solids of revolution; Green's theorem for planar areas.
- Line/circle and circle/circle intersection, point-in-polygon by ray crossing (computational geometry folklore, e.g.
  O'Rourke, *Computational Geometry in C*).
