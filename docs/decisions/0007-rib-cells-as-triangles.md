# 0007. Rib and tip-closure cells are triangles

**Decision.** From 0.9.3, rib faces and the tip closure are meshed as triangles by
an unstructured fill with boundary recovery. Their boundary nodes are the quad
sheet's nodes. This supersedes the 0.9.2 rule of one polygon of at least five
nodes per cell.

**Why.** A consumer that split polygons itself could put a triangle outside its
cell at a reflex corner (G19). Triangles from the kernel are valid by
construction, and `tests/test_rib_triangles_cover_cells.py` checks they neither
fold nor overlap.
