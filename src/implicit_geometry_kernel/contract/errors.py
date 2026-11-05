"""Shared error types for the section geometry/mesh contract.

Lives on the contract seam rather than inside either package so the mesh worker
can raise and the FE/sizing side can catch without importing geometry internals.
The geometry, mesh and FE kernels stay separate workers; this is the vocabulary
they share for "the requested section cannot be built".
"""

from __future__ import annotations


class InfeasibleSectionGeometry(ValueError):
    """The requested section cannot be realised at these thicknesses.

    Raised when the geometry is over-constrained rather than mis-specified: a
    cap deeper than the local cavity, a skin that closes the section, a web with
    no inner-skin intercept. The *request* is well-formed; the section just does
    not exist.

    Consumers that search over thicknesses (the station sizing loop) should
    catch this and score a failed trial so the search can back off, rather than
    letting it kill the worker. Subclasses ``ValueError`` so pre-existing
    ``except ValueError`` handlers keep working.
    """


__all__ = ["InfeasibleSectionGeometry"]
