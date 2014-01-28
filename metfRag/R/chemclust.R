####################################################
##
## A set of functions to cluster sets of compounds 
## by chemical similarity, and optionally to
## calculate the Maximum common substructure (MCSS)
## and on top of everything an interactive plot
##
## Input can be 1) MetFrag result SDFs, or
## DataFrames with SMILEs/InChIs and optionally scores.
##

#' Chemical Clustering of MetFrag and MetFusion results SDF files
#' 
#' Cluster a set of compounds ranked with MetFrag and MetFusion
#' 
#' @aliases plotMetClust
#' @usage plotMetClust(mols, sdf, hclust, scoreprop="Score", idprop="DatabaseID")
#' @param hclust An object of class hclust which describes the tree produced by a prior clustering process. 
#' @param mols The a list of rCDK \code{mols}. 
#' @param smiles A character vector of SMILES codes. 
#' @param filename A filename to an SD File. If boths \code{mols} 
#' and \code{filename} are given, filename is ignored. Files are read by \code{cdk::load.molecules()}
#' @param type If a \code{filename} was given,  the \code{type} will be passed to \code{cdk::load.molecules()}. 
#' The default value \code{auto} uses the filename extension.
#' @param scoreprop The name of the property of the molecules where the score is kept
#' @param idprop The name of the property of the molecules where the database ID is kept
#' @import squash
#' @importFrom tools file_ext
#' @examples 
#'        smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
#'                    'c1ccccc1Cc1ccccc1',
#'                    'C1CCC1CC(CN(C)(C))CC(=O)CC')
#'        mols <- parse.smiles(smiles)
#'        dummy <- mapply(set.property, mols, "Score", c(1,2,3,4,5))
#'        dummy <- mapply(set.property, mols, "DatabaseID", c("C1", "C2", "C3", "C4", "C5"))
#'        plotMetClust(mols)
#' 
#' @export

plotMetClust <- function(mols=NULL, filename=NULL, type=c("auto", "sdf", "smi"), 
                         hclust=NULL, scoreprop="Score", idprop="DatabaseID")
{
  
}

hclust.mols <- function(mols=NULL, smiles=NULL, filename=NULL,
                        scoreprop="Score", idprop="DatabaseID") 
{

  if (!missing(mols)) {
    ## Placeholder
    ## Nothing tbd.
  } else if (!missing(smiles)) {
    mols <- parse.smiles(smiles)    
  } else if (!missing(filename)) {
    mols <- load.molecules(filename)
  }  else {
    stop("no compounds provided")
  }

  ## Prepare compounds
  invisible(sapply(mols, do.typing))
  invisible(sapply(mols, do.aromaticity))    

  
  ## Calculate fingerprints and calculate dists
  fps <- lapply(mols, get.fingerprint, type = "extended")
  fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
  fp.dist <- 1 - fp.sim
  cluster  <- hclust(as.dist(fp.dist))

  # Put additional (chemical) information into hclust object
  cluster$mols <- mols
  cluster$labels <- sapply(mols, get.title)
  cluster$scores <- getScores(mols, scoreprop)
  cluster$ids <- getDatabaseIDs(mols, idprop)
  
  return(cluster)
}  

# h1 <- hclust.mols(filename="/vol/R/rguha/cdkr/data/BL.sdf")
# h2 <- hclust.mols(mols=NULL, filename=NULL, scoreprop="Score", idprop="DatabaseID")
# h3 <- hclust.mols(mols=NULL, filename=NULL, scoreprop="Score", idprop="DatabaseID")

getScores <- function(x, scoreprop="Score"){
  return (max(sapply(x, function(x) as.numeric(get.property(x, scoreprop)))))
}

getDatabaseIDs <- function(x, idprop="DatabaseID") {
  return (sapply(1:length(x), function(i) {
    id <- get.property(x[[i]], idprop)
    if (is.na(id) || id == 'NA')
      return (i)
    else
      return (get.property(x[[i]], idprop))
  }))
}


#' Chemical Clustering
#' 
#' Cluster a set of compounds and plot the results
#' 
#' @aliases plotCluster
#' @usage plotCluster(mols, scoreprop="Score", idprop="DatabaseID")
#' @param mols The a list of rCDK \code{mols}
#' @param scoreprop The name of the property of the molecules where the score is kept
#' @param idprop The name of the property of the molecules where the database ID is kept
#' @param ... the remaining parameters are passed down to \code{dendromat()}
#' @import squash
#' @examples 
#'        smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
#'                    'c1ccccc1Cc1ccccc1',
#'                    'C1CCC1CC(CN(C)(C))CC(=O)CC')
#'        mols <- parse.smiles(smiles)
#'        dummy <- mapply(set.property, mols, "Score", c(1,2,3,4,5))
#'        dummy <- mapply(set.property, mols, "DatabaseID", c("C1", "C2", "C3", "C4", "C5"))
#'        plotCluster(mols)
#' 
#' @export

plotCluster <- function(mols, scoreprop="Score", idprop="DatabaseID", ...) {
  
  clusters  <- list()
  
  ## Filter most promising candidates
  scores <- as.numeric(sapply(mols, get.property, scoreprop))
  mols <- unlist(mols)[scores>=0.1]
  score.explained <- scores[scores>=0.1]
  
  ## Bail out if only 1 mol left
  if (length(mols) <2) {
    warning("No mols left after filtering")
    plot(1)
    return(NULL)
  }
  
  ## Prepare compounds
  i <- 1
  invisible(sapply(mols, do.typing))
  invisible(sapply(mols, do.aromaticity))
  
  ##
  fps <- lapply(mols, get.fingerprint, type = "extended")
  fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
  fp.dist <- 1 - fp.sim
  cluster  <- hclust(as.dist(fp.dist))

  ## prepare for plotting
  classlabel   <- as.factor(paste(getDatabaseIDs(mols, idprop),
                                  cutree(cluster, h=0.2), sep=" "))
  clusterscore <- tapply(mols, classlabel, getScores)
  score.cmap <- makecmap(score.explained,
                         n = 8,
                         colFn = colorRampPalette(c("black","red")))
  clusterNumber <- length(levels(classlabel))
  cluster.cmap <- makecmap(as.numeric(classlabel),
                           n = clusterNumber)
  cluster.mat <- data.frame(Score=cmap(score.explained, score.cmap))
  
  
  dendromat(cluster, cluster.mat,
            ylab = 'Tanimoto distance', ...)
  
  legend(x=0.92,y=0.7,
         legend=score.cmap$breaks[-1], # exclude.lowest
         fill=score.cmap$colors,
         title="Score >",
         xpd=NA, bty = "n")
  
}


mcss <- function(x) {
  l <- length(x)
  
  if (l==1)
    return (x[[1]])
  
  consensus <- x[[1]]  
  for (i in 2:l) {
    consensus <- get.mcs(consensus, x[[i]])
  }
  return (consensus)  
}
