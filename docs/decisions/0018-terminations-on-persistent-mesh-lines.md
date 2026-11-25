# 0018. Terminations ride persistent mesh lines; count-changing bands are an opt-in later step

**Decision.** In the partition mesher (0015), a feature that starts or stops
along the span (a web, a cap edge that moves to the next web, a partition
boundary, a sub-partition band edge) does not change the mesh's node count. Its
mesh line exists at every station of the blade:

* where the feature is, the line is on it;
* where it is not, the line is held between the same two neighbouring lines, at
  the same fraction of the gap, as at the nearest station that has the feature;
* a partition boundary that sits on a web foot rides that foot's line, so a cap
  edge that jumps from one web to the next moves between two existing lines;
* partition boundaries are followed along the span by the kinds of partition on
  either side and by the web-foot line they follow, not by name (names
  renumber as webs and caps come and go); a cap edge that crosses a web foot
  starts a new line on the far side instead of crossing one;
* a web that ends is a plain all-quad strip that stops on its end station, taken
  from the section just inside its span;
* each skin quad takes its partition and properties from the lower station of
  its bay.

The mesher refuses (and `method="auto"` falls back to the ruled-strip mesh)
where lines would cross or meet at some stations only. Count-changing
transition bands (2:1 or 3:1 quad templates that retire a held line inside one
bay) are not the default. They may come later as an opt-in coarsening, applied
only where a held line makes elements thinner than a quality limit, and only
after tests show them against this mesh, which stays available unchanged.
Maintainer's request 2026-11-25: "whatever is best for the long term
generality".

**Why.** One rule covers any number of terminations in any order, on any run,
with no special cases per feature kind. The skin stays all-quad with every node
regular (K7, 0005), node rows keep the same identity from root to tip so
consumers can share node ids along the span, and sub-partition edges (0016) are
just more lines. Count-changing bands bring irregular nodes, skewed transition
elements, a template per ratio and per neighbour, and node rows whose identity
changes mid-span; that is more code and more ways to fail for a gain (fewer
elements far from a termination) that only matters on blades with many
terminations. The cost of persistent lines is element count and, where a held
line squeezes against a neighbour, thin elements (aspect up to about 60 on
`5D-F-LE-0R` at h = 0.25 with its 5 to 7 mm aft web); the opt-in coarsening is
the answer to that, not a different default.

**Follows.** Replaces the "structured transition bands" line in 0015's
"Follows" with this rule. Skin partitions that appear, vanish or renumber along
the span use the same held-line rule. Newton junction nodes (0017) move feet on their
lines, not the line count.
