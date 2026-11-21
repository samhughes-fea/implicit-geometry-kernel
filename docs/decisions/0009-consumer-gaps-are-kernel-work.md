# 0009. Geometry a consumer computes on top of the kernel is kernel work

**Decision.** A consumer hands all of its geometry to the kernel. Anything an
adapter still computes, converts, approximates or refuses is a kernel to-do,
listed in the consumer adapter gaps table of
[the plan doc](../internal/REFACTOR_PLAN.md) until the adapter can drop it.

**Why.** Otherwise geometry logic drifts back into consumers and two tools give two
answers for the same blade.

**Follows.** The fix is generic: it never names or special-cases the consumer
that found the gap.
