# 0004. Cap and skin share the skin mid-surface; eccentricity is applied downstream

**Decision.** Shell elements over a spar cap sit on the skin mid-surface, the same
surface as the neighbouring skin. The cap's offset from that surface (half the cap
thickness) is carried as data on each element and applied by consumers in their
constitutive law. Caps are not separate shells tied to the skin.

**Why.** Tied separate shells add constraints and complexity across every reference
blade for very little accuracy (merged and tied section stiffness agree within
about 0.1 % in EA and EI, and 1.8 % in GJ).

**Follows.** Kernel nodes never step at cap edges. The shell graph keeps the cap
partitions so laminates can be bound to them.
