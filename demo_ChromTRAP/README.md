# Demo ChromTRAP analysis

This folder contains a minimal demo for the ChromTRAP analysis using a downsampled Seurat object.

## Files

- `SO_demo.rds`: downsampled Seurat object containing BLA BA and LA excitatory neurons
- `demo_ChromTRAP.Rmd`: executable R Markdown source
- `demo_ChromTRAP.md`: rendered GitHub-friendly report

## Run

```r
rmarkdown::render("demo_ChromTRAP.Rmd")