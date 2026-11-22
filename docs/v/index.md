# V · Validation

The V manual is the evidence that the kernel computes what the
[R manual](../r/index.md) derives. Each V document evaluates the code against a
closed form, an independent computation or a published reference, states the
tolerance, and names the test that keeps the check running in CI. Each R document
links to the V document that checks it, and each V document links back.

## Documents

| Key | Title | Checks | Evidence today |
|---|---|---|---|
| V1 | Analytic sections | Circle, ellipse, thin rectangular box and flat plate against closed-form area, centroid and second moments; frame transforms against direct integration [R1] | [frames example 01](../../examples/frames/01_frames_and_twist/README.md), [second-moment rotation test](../../tests/frames/unit/test_second_moments_rotate_with_points.py) |
| V2 | Integration convergence | Quadtree integrals against the shoelace formula as depth increases [R7] | [integration example 01](../../examples/integration/01_quadtree_vs_shoelace/README.md) |
| V3 | Mesh convergence and element quality | Aspect ratio, skew and conformity of the shell mesh under refinement [R9] | [shell mesh example](../../examples/blade/15_shell_mesh/README.md), [mesh health tests](../../tests/test_mesh_health.py) |
| V4 | Golden numbers | What each snapshot pins and why it may move | [golden snapshots](../../tests/golden/), `python scripts/snapshots.py` |
| V5 | Reference blades | Summary figures and properties of public reference blades | [reference blades example](../../examples/blade/08_reference_blades/README.md) |
| V6 | Consumer compatibility | Downstream packages built and tested against each change | `consumer` job in [tests.yml](../../.github/workflows/tests.yml) |

All six documents are planned. Until they are written, the evidence column lists
the tests and examples that already make each check.
