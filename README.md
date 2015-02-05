MetFragR
========

R package for MetFrag

Installation
------------

##### Requirements

- packages built was successfully tested with R version 3.0.1 and rcdk_3.2.8.3 <br>
- first install KEGGREST dependency from Bioconductor (in R): <br>
```R
> source("http://bioconductor.org/biocLite.R")
> biocLite("KEGGREST")
```

##### Installing the package locally
- checkout the MetFragR git repository
- build the package (on command line): <br>
```bash
# R CMD check metfRag
# R CMD build metfRag
```
- install the package (in R): <br>
```R
> install.packages("metfRag",repos=NULL,type="source")
> library(metfRag)
```

##### Installing the package directly from github
```R
> library(devtools)
> install_github("c-ruttkies/MetFragR/metfRag")
```

If you get an error like "ERROR: loading failed for 'i386'" on windows
due to rJava issues, disable multiarch builds:
```R
install_github("c-ruttkies/MetFragR/metfRag", args="--no-multiarch")
```

