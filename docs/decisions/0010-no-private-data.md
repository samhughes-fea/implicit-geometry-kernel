# 0010. No confidential blade data in this repository

**Decision.** This repository holds no inputs for non-public blades: no aerofoil
coordinates, layups, planforms or material cards. Reference blades are read from a
separate catalogue checkout, found through `REFERENCE_BLADES_DIR` or a git-ignored
`reference_blades/` folder, and tests and examples skip without it. Private
projects are not named anywhere in the repository.

**Why.** Some reference blades are confidential, and the repository is public.

**Follows.** Figures and summary results of public reference blades may be
committed. `tests/test_consumer_names.py` checks for private names.
