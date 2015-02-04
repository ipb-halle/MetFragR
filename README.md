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
- checkout the MetFragR git repository <br>
command line:
```bash
R CMD check metfRag
R CMD build metfRag
```
R:
```R
> install.packages("metfRag",repos=NULL,type="source")
> library(metfRag)
```

##### Installing the package directly from github
R:
```R
library(devtools) <br>
install_github("c-ruttkies/MetFragR/metfRag")
```
