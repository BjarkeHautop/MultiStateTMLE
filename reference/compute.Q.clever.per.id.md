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
# S = 1, single terminal process "D": Q_t collapses to the discrete-time
# cumulative incidence 1 - prod(1 - h_s).
h <- c(0.1, 0.2, 0.3)
dt_id <- data.table::data.table(
  state = c(1L, 1L, 1L),
  state.row.index = c(1L, 1L, 1L),
  P.D.1 = h
)
states <- data.table::data.table(state = 1L)

compute.Q.clever.per.id(
  dt_id = dt_id,
  states = states,
  process.types = list(D = "terminal")
)
#>    state state.row.index P.D.1     Q clever.Q.D0 clever.Q.D1 clever.Q.D   Q.1
#>    <int>           <int> <num> <num>       <num>       <num>      <num> <num>
#> 1:     1               1   0.1 0.496        0.44           1       0.56 0.496
#> 2:     1               1   0.2 0.440        0.30           1       0.70 0.440
#> 3:     1               1   0.3 0.300        0.30           1       0.70 0.300
```
