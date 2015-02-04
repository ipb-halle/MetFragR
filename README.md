MetFragR
========

R package for MetFrag

Installation
------------

##### Requirements

- packages built was successfully tested with R version 3.0.1 and rcdk_3.2.8.3 <br>
- first install KEGGREST dependency from Bioconductor: <br>
R:
```R
> source("http://bioconductor.org/biocLite.R")
> biocLite("KEGGREST")
```

##### Installing the package locally
- checkout the MetFragR git repository
- build the package: <br>
command line:
```bash
R CMD check metfRag
R CMD build metfRag
```
- install the package: <br>
R:
```R
> install.packages("metfRag",repos=NULL,type="source")
> library(metfRag)
```

##### Installing the package directly from github
R:
```R
library(devtools)
install_github("c-ruttkies/MetFragR/metfRag")
```
