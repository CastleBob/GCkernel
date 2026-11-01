# Independent IGES manifold-solid intake

`iges::read` and `read_file` own IGES 5.3 ASCII records, directory/parameter
references, units and explicit B-rep construction. Independently written from
the US PRO specification, sections 2.2.3/2.2.4 and entity tables; no translator
dependency or third-party implementation source is used.

The first domain is 1..1024 type-186 solids with one closed type-514 shell,
type-510 faces, type-508 oriented loops, shared type-504 edges and type-502
vertices. Identity comes from directory/index pairs, never coordinate sewing.
Geometry: lines (110), arcs (100), points (116), directions (123), proper rigid
transforms (124), parameterized planes/cylinders (190/192), single clamped
degree-1..3 positive-weight rational curves/surfaces (126/128). Exact homogeneous
elevation and parameter restriction preserve source geometry. IGES controls are
U-fast; neutral patch controls are U-major.

Explicit spline pcurves must cover every edge use, one bounded line/cubic per
edge, normalized from the original surface domain and certified against the
physical edge. Contradictions fail. Periodic/closed/multispan IGES surfaces,
analytic pcurves, composite curves, voids, loose sheets, construction/drawing
entities, unreferenced entities and unknown types fail explicitly. No geometry
is omitted or replaced. This is bounded initial coverage, not general CAD parity.

The parser validates 80-column records, section order/sequences and termination
counts, directory ownership/ranges, Hollerith strings, delimiters, finite reals
(including D exponents), and complete consumption. Limits: 64 MiB, 500000 entities,
four million parameters per entity. Unit flags/names and model-to-real scale
determine inches; only documented format defaults are applied. File resolution
does not increase admission tolerance. Strict neutral, topology, connectedness,
winding and positive-mass checks precede any return. No file writes or partial results.

Op 22 request: `[1,pathByteCount,UTF8 bytes]`. Reply: `[1,count,(directoryId,
nameByteCount,UTF8 name,sourceToInches,volume,bodyLength,body)...]`.
`KernelIges.Read` is structural transport only. Production intake remains EyeShot.
Analytic fixtures cover mass, units/model scale, rational geometry and placement;
native comparisons are separate evidence. Malformed/unsupported input must fail.

Specification: https://shocksolution.files.wordpress.com/2017/12/iges5-3_fordownload.pdf

Independent-intake caveat: the model-scale .5 millimetre spline box is 24 cubic
inches by the authored dimensions and IGES 5.3 section 2.2.4.3.13. Existing native
intake returns 3 cubic inches, leaving model scale unapplied. This disagreement
is characterized separately; it is not counted as native agreement. The same
rational geometry with model scale 1 is a separate comparison fixture.
