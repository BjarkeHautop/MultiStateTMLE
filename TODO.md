# Performance:

- R/prepare.initial.R: repeated `grep()`/`setdiff()`/`unique()` recomputation on the same vectors inside inner loops. Needs a closer look with tests before touching.
