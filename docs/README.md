# Reproducible report

`reproducible_report.Rmd` recomputes every quantity in the commentary in order,
with an explanation of what each step does and why. It reproduces the original
Article's headline results before departing from them.

It duplicates the analysis in `R/` rather than replacing it: the numbered
scripts are the deposit, and this document is a narrated verification of them.

Render with:

```r
rmarkdown::render("docs/reproducible_report.Rmd")
```

This requires the Code Ocean capsule files in `data/raw/` (see
`data/raw/README.md`). Commit the rendered `reproducible_report.html` alongside
the source so readers can read it without running anything.
