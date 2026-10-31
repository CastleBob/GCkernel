# STL intake contract

Independently written ASCII/binary STL parsing owns file syntax, explicit unit
conversion and detached facets. `read`/`read_file` require a finite positive
source-length-to-inch factor; the format contains no trustworthy unit declaration.
There is no guessed scale. Binary classification uses the exact 84+50*n byte
length before parsing, including headers starting with `solid`. ASCII requires
one complete solid block. Nonfinite/degenerate facets, malformed/trailing records,
nonzero binary attributes and inconsistent nonzero normals are explicit errors.
A zero declared normal is explicitly unspecified; vertex winding is authoritative.
No input normal or vertex is modified to make a failed solid succeed.

`to_solid` is an explicit separate operation. Vertices share topology only when
coordinates are exactly equal (signed zeros compare equally); there is no tolerance
welding. Every edge must occur once each way, each vertex link must be one cycle,
and the shell must be connected and outward oriented. A bounding-box sweep and
separating-axis certificates reject crossing/overlapping facets or unresolved
contact. Shared topology allows only that exact vertex/edge contact. The result
retains planar triangular faces; no analytic surface recognition or smoothing.
Open/disconnected meshes remain available through mesh intake, never substituted
for a failed solid request. Multiple-shell/cavity solid conversion remains open.

Limits: 64 MiB input, 100000 facets, positive finite facet double area above 1e-14
square inches, two million intersection candidate pairs. Inconclusive floating
point geometry is refused. Stored nonzero normals must be unit and agree with
winding within 1e-4 (binary values are float32); this is not a position tolerance.
There is no third-party code or dependency.

Worker operation 18 request:
`[1, mode (0 mesh / 1 solid), sourceLengthToInches, pathByteCount, UTF8 bytes]`.
Mesh reply: `[1,0,facetCount,(declaredNormalXYZ,Axyz,Bxyz,Cxyz)...]`.
Solid reply: `[1,1,volume,bodyValueCount,neutralBody...]`.
C# KernelStl transports detached data only. File/document writes and production
import UI are outside this owner. ReadSolid never retries with ReadMesh.

Primary public format descriptions reviewed 2026-10-31:
- https://www.loc.gov/preservation/digital/formats/fdd/fdd000505.shtml
- https://www.loc.gov/preservation/digital/formats/fdd/fdd000506.shtml

Regression boundary: independently authored tetrahedra in both encodings, analytic
volume at inch/mm scales, shared coplanar triangles, large translations, STEP
export, truncation/nonfinite/color refusal, open/reversed/cracked/disconnected
meshes, pinched vertices, overlapping and intersecting facets, and transport.
