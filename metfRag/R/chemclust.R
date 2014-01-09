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



getscores <- function(x, scoreprop="Score"){
  return (max(sapply(x, function(x) as.numeric(get.property(x, scoreprop)))))
}

getDatabaseIDs <- function(x, idprop="DatabaseID") {
  return (sapply(1:length(x), function(id) {
    if (get.property(x[[id]], idprop) == 'NA')
      return (id)
    else
      return (get.property(x[[id]], idprop))
  }))
}


#' Chemical Clustering
#' 
#' Cluster a set of compounds and plot the results
#' 
#' @aliases plotCluster
#' @usage plotCluster(mols, scoreprop="Score", idprop="DatabaseID")
#' @param w The a list of rCDK \code{mols}
#' @param scoreprop The name of the property of the molecules where the score is kept
#' @param idprop The name of the property of the molecules where the database ID is kept
#' @examples 
#'        smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
#'                    'c1ccccc1Cc1ccccc1',
#'                    'C1CCC1CC(CN(C)(C))CC(=O)CC')
#'        mols <- parse.smiles(smiles)
#'              plotCluster(mols)
#' 
#' @export

plotCluster <- function(mols, scoreprop="Score", idprop="DatabaseID") {
  
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
  classlabel   <- as.factor(paste(getDatabaseIDs(mols, idprop),
                                  cutree(cluster, h=0.2), sep=" "))
  clusterscore <- tapply(mols, classlabel, getscores)
  score.cmap <- makecmap(score.explained,
                         n = 8,
                         colFn = colorRampPalette(c("black","red")))
  clusterNumber <- length(levels(classlabel))
  cluster.cmap <- makecmap(as.numeric(classlabel),
                           n = clusterNumber)
  cluster.mat <- data.frame(Score=cmap(score.explained, score.cmap))
  
  
  dendromat(cluster, cluster.mat,
            ylab = 'Tanimoto distance', main = "")
  
  legend(x=0.96,y=0.7,
         legend=c("0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95"),
         fill=score.cmap$colors,
         title="Score >",
         xpd=NA, bty = "n")
  
}
