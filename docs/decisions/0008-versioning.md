# 0008. Releases stay on 0.9.x; tags are the maintainer's

**Decision.** Every release is `0.9.N`. There is no 0.10.0, and 1.0.0 is cut only
when the maintainer decides the public API is settled. Only the maintainer creates
tags and GitHub releases.

**Why.** The public API is still being worked out; staying pre-1.0 says so, and
patch numbers keep consumers' pins simple.

**Follows.** Pull requests bump `__version__` and the changelog heading;
`release.yml` refuses a tag that does not match. Removals follow
[API stability](../conventions/api/STABILITY.md).
