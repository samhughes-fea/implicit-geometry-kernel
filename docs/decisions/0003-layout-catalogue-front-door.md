# 0003. The layout catalogue is the front door; ribs are a spanwise axis

**Decision.** The parametric layout catalogue (`WX-Y-Z`, 40 layouts) is how most
callers describe a section. Each id compiles to a list of member contributions.
Ribs, bulkheads and stringers are a second, spanwise axis of a concept (`NR` and
`SpanwiseScheme`), not another letter in the code. The public name for webs is
`ShearWeb`.

**Why.** A short code that names a real structural concept, with the spanwise
scheme kept separate so the two vary independently.
