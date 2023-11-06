# RNAseq example
Showcasing bulk RNAseq processing and analysis. This analysis is almost completely containerized (see Dockerfile) and reproducible (data are downloaded straight from SRA). The exception is that some of the R packages (e.g., DESeq2) were difficult to conerce into the image. Therefore, the analysis and figures, done in R, were completed in RStudio in an .Rmd file.
