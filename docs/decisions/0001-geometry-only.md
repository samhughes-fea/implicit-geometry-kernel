# 0001. The kernel is geometry only

**Decision.** The kernel describes blade and wing section geometry and its span
loft. Materials, plies, stiffness, loads and solvers live in downstream packages.
There is no separate "blade model" package between them (step 9 was declined).

**Why.** One geometry source for every downstream tool, and a package small enough
to keep exact and well tested.

**Follows.** Thickness is a law in metres, not a layup. Seams that let downstream
tools attach materials (partitions, named edges, forced stations) live here; the
materials do not. See [SCOPE.md](../../SCOPE.md).
