# Compute Q and clever covariates per id using a discrete-state backward recursion

The function expects a time-ordered `dt_id` (rows for a single id) that
contains hazard columns named like `P.<name>.<s>` for each discovered
process `<name>` and state index `s` in `1..S`. It uses the
product-space `states` table to map state indices and to build successor
mappings.

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

  data.table with time-ordered rows for a single id. Must contain
  columns for hazards of the form P..~~.~~

- states:

  data.table enumerating the product state space. Must have S rows and
  columns for each stateful process (named by process) plus a column
  `state` giving indices 1..S.

- process.types:

  named character vector or list describing each discovered process
  type. Allowed values include: "terminal", "recurrent", "one.jump",
  "state-with-atrisk". If omitted, a small heuristic is used.

- P.prefix:

  prefix used for hazard columns in dt_id (default "P."). The function
  looks for columns matching `^P.prefix<name>.<s>$`.

- state.idx.col:

  name of the column in dt_id that holds the state index (value in
  1..S). Default "state".

- parameter:

  character: name of the target process/outcome. If "target" (default),
  the first discovered terminal/outcome name is used.

- process.deltas:

  optional numeric vector (length = number of discovered processes). Not
  used internally for mapping, only validated if provided.

- compute.clever:

  logical; whether to compute and append clever.Q.0 / clever.Q.1
  (default TRUE).

- browse:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html) at the start of the
  backward recursion.

- browse2:

  logical; if `TRUE`, drop into
  [`browser()`](https://rdrr.io/r/base/browser.html) at a secondary
  breakpoint inside the recursion.

- get.years.lost:

  logical; whether to also compute a years-lost-type summary.

- years.lost.block.size:

  block size used when computing the years-lost summary (when
  `get.years.lost = TRUE`).

- clever.by.state:

  logical; whether to compute clever covariates by state rather than
  pooled.

## Value

data.table: `dt_id` augmented with column `Q` and columns
`clever.Q.<name>0` and `clever.Q.<name>1` for each discovered name.

## Examples

``` r
# compute.Q.clever.per.id(dt_id = some_dt_for_one_id, states = depend.matrix, process.types = process.types)
```
