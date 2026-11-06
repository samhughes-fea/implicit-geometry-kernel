# 11 · Partial-span free-end stringer

A `Stringer` is a branch that stands off one skin. `extent` limits it to part of
the span; those ends are `TopologyEvent`s and `forced_stations()`. Inside the
interval the branch appears on `shell_graph().branches`; outside it is absent.
It is not lofted with the skin loop.

```bash
python examples/11_partial_span_stringer/partial_span_stringer.py
```
