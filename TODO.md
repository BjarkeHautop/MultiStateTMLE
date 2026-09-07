# Breaking changes desired

- Use snake_case.

- Clean up arguments (and defaults).

# CRAN compliance

- Can't set number of cores > 1 by default.

# Performance:

- R/prepare.initial.R: repeated `grep()`/`setdiff()`/`unique()` recomputation on the same vectors inside inner loops. Needs a closer look with tests before touching.
