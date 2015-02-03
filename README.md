MetFragR
========

R package for MetFrag

Installation
------------

\# R CMD check metfRag <br>
\# R CMD build metfRag <br>
\# R <br>
\> # first install KEGGREST dependency <br>
\> source("http://bioconductor.org/biocLite.R") <br>
\> biocLite("KEGGREST") <br>
\> # install the r package <br>
\> install.packages("metfRag",repos=NULL,type="source") <br>
\> library(metfRag) <br>

built with R version 3.0.1 and rcdk_3.2.8.3

or directly with the github install command:<br>

\# R <br>
\> library(devtools) <br>
\> install_github("c-ruttkies/MetFragR/metfRag")
