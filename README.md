MetFragR
========

R package for MetFrag

Installation
------------

\# R CMD check metfRag <br>
\# R CMD build metfRag <br>
\# R <br>
\> # first install KEGGREST dependency
\> source("http://bioconductor.org/biocLite.R")
\> biocLite("KEGGREST")
\> # install the r package
\> install.packages("metfRag",repos=NULL,type="source") <br>
\> library(metfRag) <br>

built with R version 3.0.1 and rcdk_3.2.8.3

or directly with the github install command:

\# R <br>
\> library(devtools)
\> install_github("c-ruttkies/MetFragR/metfRag")
