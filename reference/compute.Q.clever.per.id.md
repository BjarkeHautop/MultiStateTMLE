# Compute Q and clever covariates per id (deprecated)

`compute.Q.clever.per.id()` is deprecated as of MultiStateTMLE 0.1.0 and
is now internal; it is called by
[`tmle_alpha()`](https://github.com/BjarkeHautop/MultiStateTMLE/reference/tmle_alpha.md).

## Usage

``` r
compute.Q.clever.per.id(
  dt_id,
  states,
  process.types = NULL,
  P.prefix = "P.",
  state.idx.col = "state",
  parameter = "target",
  process.deltas = NULL,
  compute.clever = TRUE,
  browse = FALSE,
  browse2 = FALSE,
  get.years.lost = FALSE,
  years.lost.block.size = 10,
  clever.by.state = FALSE
)
```

## Arguments

- dt_id:

  A `data.table` with the time-ordered rows of a single id, with hazard
  columns named `<P.prefix><name>.<s>` for each process `<name>` and
  state index `s`.

- states:

  A `data.table` enumerating the product state space, with a column per
  stateful process and the state index column.

- process.types:

  Named list of process types (`"terminal"`, `"recurrent"`,
  `"one.jump"`).

- P.prefix:

  Prefix of the hazard columns.

- state.idx.col:

  Name of the state index column.

- parameter:

  Name of the target process; `"target"` picks the first terminal
  process.

- process.deltas, browse, browse2:

  Ignored.

- compute.clever:

  Logical. Whether to compute clever covariates.

- get.years.lost, years.lost.block.size:

  Whether to compute the expected time lost, and the block size of its
  approximation.

- clever.by.state:

  Logical. Whether to compute state-specific clever covariates.

## Value

`dt_id` with `Q`, `Q.<s>` and `clever.Q.<name>` columns added.
