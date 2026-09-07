# Breaking changes desired

- Use snake_case.

- Clean up arguments (and defaults).

# CRAN compliance

- Can’t set number of cores \> 1 by default.

- Add examples and tests.

# Performance:

- R/prepare.initial.R: repeated
  [`grep()`](https://rdrr.io/r/base/grep.html)/[`setdiff()`](https://rdrr.io/r/base/sets.html)/[`unique()`](https://rdrr.io/r/base/unique.html)
  recomputation on the same vectors inside inner loops. Needs a closer
  look with tests before touching.
