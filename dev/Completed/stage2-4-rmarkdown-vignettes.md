# Stage 2.4 R Markdown vignette migration

Stage 2.4 migrates the user-facing vignette sources from Sweave-style `.Rnw` files to R Markdown `.Rmd` files.

The migration keeps the same broad examples while making the vignette sources easier to maintain, easier to diff, and more consistent with modern R package documentation workflows.

Generated vignette outputs and the old `.Rnw` sources are removed by the stage runner. The package now declares `rmarkdown` in `Suggests` because the vignette sources use `knitr::rmarkdown`.
