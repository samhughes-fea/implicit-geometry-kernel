# 0006. The tip closure is part of the skin shell, not a rib

**Decision.** Every blade closes its skin shell at z/L = 1 with a tip closure meshed
by the kernel. It carries the skin partition and material. `NR` in a layout code
counts interior ribs only: a `10R` blade has ten ribs and one tip closure.

**Why.** An open tip is not a physical blade, and treating the closure as a rib
would give it the wrong partition and laminate.
