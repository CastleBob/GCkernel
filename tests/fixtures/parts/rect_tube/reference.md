# Rect Tube.gcmprt: kernel vs an Eyeshot step-by-step reference

Reference = Eyeshot Brep.Difference of the same captured tools in order (not GCmaker's conditioned
production path). Volumes: kernel MEASURE (exact or quadrature) and the reference's mesh at chord 1e-6 x size.
Probes: 13^3 grid over the step's tool box grown .3, kernel vs reference.

| step | tool | kernel | kernel volume | Eyeshot reference volume (mesh) | probes differing |
|---|---|---|---|---|---|
| 1 | Left end condition operation 1 tool 1 | complete | 61.804389648 (closed form) | 61.804262796 | 0 of 2197 |
| 2 | Hole_1 | complete | 58.662796995 (closed form) | 58.662702784 | 0 of 2197 |
| 3 | Extrude1 (Cut) | complete | 55.725573602 (closed form) | 55.725486031 | 0 of 2197 |
| 4 | Extrude2 (Cut) | complete | 54.783787134 (closed form) | 54.783702062 | 0 of 2197 |
| 5 | Extrude3 (Cut) | INCOMPLETE: stopped at step 5 of 5 (Extrude3 (Cut)): not supported: no cardinal section chart for the planar cutter and stock cylinders Showing the body after step 4. Base imported from GCmaker. | - | 54.342254292 | - |

GCmaker's own final solid (all operations, production path): mesh volume 54.342267306.
