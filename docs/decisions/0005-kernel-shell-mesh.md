# 0005. The kernel ships the analysis shell mesh; consumers do not remesh

**Decision.** `shell_mesh` is the production mesh. Consumers use its nodes and
elements as they are. Gmsh, where a consumer keeps it, is for comparison and
refinement studies only. The kernel has no Gmsh driver.

**Why.** One mesh that conforms to the kernel's own partitions, so every consumer
analyses the same discretised geometry.

**Follows.** Skin and web mid-surfaces are all-quad mapped strips between span
stations. Refinement controls (edge size, leading-edge and spanwise refinement)
live in the kernel.
