# kfltCNV
A tool for calling CNV (Copy number variations) by the Kalman Filter

## Install
It can be installed from github:

```
devtools::install_github('sanadamakomi/kfltCNV')
```

If you wants to build the vignettes, run:

```
devtools::install_github('sanadamakomi/kfltCNV', 
  build = TRUE, build_opts = c("--no-resave-data"), force = TRUE)
```

## Test 

There are test Data in the directory 'kfltCNV/inst/extdata', run scripts in R console:

```
library(kfltCNV)
testBam <- system.file("extdata", 'testSample.bam', package = "kfltCNV")
controlBam1 <- system.file("extdata", 'controlSample.bam', package = "kfltCNV")
controlBam2 <- system.file("extdata", 'controlSample2.bam', package = "kfltCNV")
controlBam <- c(controlBam1, controlBam2)
bed <- system.file("extdata", 'chr10_exome.bed', package = "kfltCNV")
annotate <- system.file("extdata", 'hg19_refGene_chr10.txt', package = "kfltCNV")
```

Pipeline to call CNVs:

```
kfltBatch(testBam, controlBamFile = controlBam, bedFile = bed)

# Annotate 
kfltBatch(testBam, controlBamFile = controlBam, bedFile = bed, annote.database = annotate, outDir = '.')
```

See details in `kfltCNV-vignette`.
