# STR-1.gcmprt: kernel vs an Eyeshot step-by-step reference

Reference = Eyeshot Brep.Difference of the same captured tools in order (not GCmaker's conditioned
production path). Volumes: kernel MEASURE (exact or quadrature) and the reference's mesh at chord 1e-6 x size.
Probes: 13^3 grid over the step's tool box grown .3, kernel vs reference.

| step | tool | kernel | kernel volume | Eyeshot reference volume (mesh) | probes differing |
|---|---|---|---|---|---|
| 1 | Left end condition operation 1 tool 1 | complete | 112.845372603 (closed form) | 112.845646699 | 0 of 2197 |
| 2 | Left end condition operation 2 tool 1 | INCOMPLETE: stopped at step 2 of 2 (Left end condition operation 2 tool 1): not supported: transverse graph coincidence is unresolved Showing the body after step 1. Base imported from GCmaker. | - | 110.414764497 | - |
| 3 | Left end condition operation 3 tool 1 | INCOMPLETE: stopped at step 2 of 3 (Left end condition operation 2 tool 1): not supported: transverse graph coincidence is unresolved Showing the body after step 1. Base imported from GCmaker. | - | 106.701206509 | - |
| 4 | Right end condition operation 1 tool 1 | INCOMPLETE: stopped at step 2 of 4 (Left end condition operation 2 tool 1): not supported: transverse graph coincidence is unresolved Showing the body after step 1. Base imported from GCmaker. | - | 101.28852158 | - |
| 5 | Right end condition operation 2 tool 1 | INCOMPLETE: stopped at step 2 of 5 (Left end condition operation 2 tool 1): not supported: transverse graph coincidence is unresolved Showing the body after step 1. Base imported from GCmaker. | - | 101.231382506 | - |
| 6 | Right end condition operation 3 tool 1 | INCOMPLETE: stopped at step 2 of 6 (Left end condition operation 2 tool 1): not supported: transverse graph coincidence is unresolved Showing the body after step 1. Base imported from GCmaker. | - | 98.800405878 | - |

GCmaker's own final solid (all operations, production path): mesh volume 98.800405878.
