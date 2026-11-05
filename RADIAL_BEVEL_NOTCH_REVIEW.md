# Radial notch across a near bevel — Codex, 2026-11-05

Branch codex/radial-bevel-notch extends b88daff within the Boolean boundary.
Implementation and contract are in src/boolean/radial_notch.rs and README.md.
The complete admission runs before domain selection. A simple outer bevel
at the crossed end adds two exact ConePair edges connected to the central
CylinderPair; the body is spliced and the old cone seam slot reused. Opposite
end topology remains byte-for-byte unchanged at the same indices.

Five new tests cover three bevel/cutter proportions, independent disk-strip
volume, neutral round trip, outward closed meshes, rotations/opposite ends,
one exact translated frame, scales 0.25/4, explicit history and six refusal
cases preserving both operands. Existing flat-notch (7), annular-window (7),
cone-bore (5) and radial-bore (4) checks also pass: 28 distinct focused tests.
The previous flat-notch refusal for a now-supported crossing was changed to
a smaller cutter confined to the cone; the admitted crossing is explicitly
tested in the new suite, not silently dropped.

For the original Pipe Coped Chamfered dimensions, the independent reference
is 5.17435183987158620 in3; kernel mass agrees within 1e-9. Actual saved-recipe
and host checks are next, not yet claimed. Source no longer requires that
bevel crossing to refuse, but repeated/coincident cuts remain unsupported.
Exact relations may refuse rounded arbitrary transforms. No guessed recovery,
new dependency, third-party implementation, native-parity waiver or pin
promotion. Claude's active checkout, GCUI and pending policy choices unchanged.

Logs: TEMP/GCKernel-radial-bevel-{first,regression,radial-regression}.log.
An initial test-file extraction error failed compilation, corrected before
running the completed suite; no runtime workaround was added.
