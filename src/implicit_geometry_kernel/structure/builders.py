"""Extension point: a builder turns design parameters into contributions.

Add a member type by subclassing :class:`MemberBuilder`; nothing else changes. Builders
run in registration order and see the contributions of the builders before them, so a
builder that joins other members (a web meeting the caps) comes after them.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass, field

from implicit_geometry_kernel.sdf.field import FeatureTable, SDFNode
from implicit_geometry_kernel.structure.contributions import Contribution


@dataclass(frozen=True)
class BuildContext:
    """What every builder may use. Builders change nothing here except ``table``.

    ``oml`` is the outer mould line (negative inside) and ``span`` the interval of
    ``z`` the body occupies. Callers subclass it to carry their own parameters.
    """

    oml: SDFNode
    span: tuple[float, float]
    table: FeatureTable = field(default_factory=FeatureTable)

    def __post_init__(self) -> None:
        if not self.span[1] >= self.span[0]:
            raise ValueError(f"span must be (z_start, z_end) with z_end >= z_start, got {self.span}")


class MemberBuilder(ABC):
    """Builds one kind of member."""

    @abstractmethod
    def build(self, ctx: BuildContext, prior: Sequence[Contribution]) -> list[Contribution]:
        """Return this builder's contributions; ``prior`` holds earlier builders' output."""


@dataclass
class BuilderRegistry:
    """Builders in the order they run."""

    builders: list[MemberBuilder] = field(default_factory=list)

    def register(self, builder: MemberBuilder) -> BuilderRegistry:
        if not isinstance(builder, MemberBuilder):
            raise TypeError(f"{builder!r} is not a MemberBuilder")
        self.builders.append(builder)
        return self

    def run(self, ctx: BuildContext) -> list[Contribution]:
        out: list[Contribution] = []
        for builder in self.builders:
            made = builder.build(ctx, tuple(out))
            out.extend(made)
        return out


__all__ = ["BuildContext", "BuilderRegistry", "MemberBuilder"]
