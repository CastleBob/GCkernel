# Opposite normalized pipe ends

Codex, 2026-11-04 night. Isolated candidate from e41761f; this commit changes
only boolean/normal_end.rs, its tests and this note. GCcore d67f9d0 is the
original red saved-recipe reproducer against174d77f.

The existing selected NormalEnd domain now admits a second opposite end when
the stock has two complete cylindrical walls, one plain planar cap and one
complete normalized conoid cap. Its radial rim recipes must match the retained
conoid exactly; its seam is one line used twice. Its axis must equal the new
axis or its negation exactly (no new angular substitution). Both caps supply
containing axial ranges using the existing pair-height bound. The new conoid
must lie strictly beyond the planar end and before the whole retained conoid
range. Global range overlap refuses even when a more refined proof might admit
a case. Holes in the input walls, recutting an existing normalized end and
missing/irregular cap topology remain outside this operation.

The retained cap surface, sense, curves and loop directions are copied exactly,
including its radial seam. The new end reuses the existing construction; no
new mass, tessellation, naming or protocol code. No retry or changed tolerance.
Input topology checks and output topology validation remain mandatory. This
is numerical geometric admission, not a formal interval-floating-point proof.

Public mathematics: [MathWorld Cylindrical Coordinates](https://mathworld.wolfram.com/CylindricalCoordinates.html),
equation17, reviewed2026-11-04. The independent radial-height integral and
separation reasoning are ours; no third-party kernel implementation/dependencies.

Verification: ten normal_end tests pass, followed by one new translated/tilted
frame test passing (eleven total). Seven were existing regressions, including
the original single-end/hole behavior. New cases cover12 two-end combinations
(twoorders, three mating-axis turns, PreserveLength-equivalent offset on/off),
exact first-cap geometry retention, volume, neutral round trip, bitwise welded
mesh closure and positive signed mesh volume. The translated frame covers both
orders. Touching/overlap/same-end recuts refuse; a radial hole after both ends
passes independent removed-volume and mesh checks. Source bodies remain intact.
Neutral Reader::direction normalizes vectors: turned-axis wire round trips are
checked within16epsilon per scalar plus topology/volume/mesh, not bit equality.

GCcore's unmodified two aggregate saved-recipe tests now pass all12 rebuilds:
eight single-end plus four two-end inch/mm and PreserveLength cases. Logs:
TEMP/GCKernel-opposite-normalized-{tests,frame}.log and
TEMP/GCcore-opposite-normalized-candidate.log.

Persistent ancestry for normalized side composites remains unavailable as before.
No saved schema, side identity, active checkout, installed pin or release change.
Claude review/adoption requested separately; no acknowledgment claimed.
