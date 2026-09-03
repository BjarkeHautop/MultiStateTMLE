# Correctness bugs

- R/make.hal.library.R:18: syntax error. The file also reference undefined globals.

- R/make.hal.R has the same shape of bug.

- R/compute.Q.clever.per.id.R:56: `which(state_processes == derived_processes)` compares vectors of possibly different length. Currently incide an `if (FALSE)` block.

- R/sim.generic.R:126,141,145: `term_deltas`/`column indexing` hard-codes "exactly one censoring process at index 0." R/sim.from.data.R:185-189 fixed this via `match(term.processes, process.order) - 1L`, so old version?

R/tmle.alpha.fun.R:173: `clever.weight` is only created if `(length(a)>0)`, but used unconditionally later in the EIC/TMLE update loop.

# Performance footguns

- R/prepare.initial.R:1385,1448,1752,1771: `depend.matrix <- rbind(depend.matrix, depend.matrix.t2)` inside for loops, growing a data.table every iteration.

- R/compute.Q.clever.per.id.R:301-308: nested loop recomputing full-length logical vectors per unique time combination, inside an outer mclapply over ids. A join/match would replace it.

- R/tmle.alpha.fun.R:379-388: double loop over process × state using `tmp.long[[col]] <- ...` (copies the column) instead of data.table's `:=` (in-place update).

- R/prepare.initial.R: repeated `grep()`/`setdiff()`/`unique()` recomputation on the same vectors inside inner loops.
