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
#' @param hclust An object of class hclust which describes the tree produced by a prior clustering process. 
#' @param mols The a list of rCDK \code{mols}. 
#' @param filename A filename to an SD File. If boths \code{mols} 
#' and \code{filename} are given, filename is ignored. Files are read by \code{cdk::load.molecules()}
#' @param type If a \code{filename} was given,  the \code{type} will be passed to \code{cdk::load.molecules()}. 
#' The default value \code{auto} uses the filename extension.
#' @param scoreprop The name of the property of the molecules where the score is kept
#' @param idprop The name of the property of the molecules where the database ID is kept
#' @import squash
#' @importFrom tools file_ext
#' @examples 
#'        library(rcdk)
#'        smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
#'                    'c1ccccc1Cc1ccccc1',
#'                    'C1CCC1CC(CN(C)(C))CC(=O)CC')
#'        mols <- parse.smiles(smiles)
#'        dummy <- mapply(set.property, mols, "Score", c(1,2,3,4,5))
#'        dummy <- mapply(set.property, mols, "Identifier", c("C1", "C2", "C3", "C4", "C5"))
#'        plotMetClust(mols)
#' 
#' @export

plotMetClust <- function(mols=NULL, filename=NULL, type=c("auto", "sdf", "smi"), 
                         hclust=NULL, scoreprop="Score", idprop="Identifier")
{
  warning("NYI")
}

#' @export
hclust.mols <- function(mols=NULL, smiles=NULL, filename=NULL,
                        scoreprop="Score", idprop="Identifier") 
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
# h2 <- hclust.mols(mols=NULL, filename=NULL, scoreprop="Score", idprop="Identifier")
# h3 <- hclust.mols(mols=NULL, filename=NULL, scoreprop="Score", idprop="Identifier")


#' Chemical Clustering
#' 
#' Cluster a set of compounds and plot the results
#' 
#' @aliases plotCluster
#' @param mols The a list of rCDK \code{mols}
#' @param scoreprop The name of the property of the molecules where the score is kept
#' @param idprop The name of the property of the molecules where the database ID is kept
#' @param k,h Scalar. Cut the dendrogram such that either exactly
#' \code{k} clusters are produced or by cutting at height \code{h}.
#' (either k or h needs to be specified)
#' @param ... the remaining parameters are passed down to \code{dendromat()}
#' @import squash
#' @examples 
#'        library(rcdk)
#'        smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
#'                    'c1ccccc1Cc1ccccc1',
#'                    'C1CCC1CC(CN(C)(C))CC(=O)CC')
#'        mols <- parse.smiles(smiles)
#'        dummy <- mapply(set.property, mols, "Score", c(1,2,3,4,5))
#'        dummy <- mapply(set.property, mols, "Identifier", c("C1", "C2", "C3", "C4", "C5"))
#'        plotCluster(mols, h=0.2)
#' 
#' @export

plotCluster <- function(mols, scoreprop="Score", idprop="Identifier", h=NULL, k=NULL, ...) {
  
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
                                  cutree(cluster, h=h, k=k), sep=" "))
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


#' Chemical Clustering
#' 
#' Calculate the Maximum Common Substructure
#' 
#' @aliases getClusterMCSS
#' @param hcluster Previously created hclust for \code{mols}
#' @param mols The a list of rCDK \code{mols}
#' @param k,h Scalar. Cut the dendrogram such that either exactly
#' \code{k} clusters are produced or by cutting at height \code{h}.
#' (either k or h needs to be specified)
#' @param which A vector selecting the clusters to be returned
#' @examples 
#'        library(rcdk)
#'        smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
#'                    'c1ccccc1Cc1ccccc1',
#'                    'C1CCC1CC(CN(C)(C))CC(=O)CC')
#'        mols <- parse.smiles(smiles)
#'        cluster <- hclust.mols(mols)
#'        clusterreps <- getClusterMCSS(cluster, mols, k=2) 
#' 
#' @export

getClusterMCSS <- function(hcluster, mols, h=NULL, k=NULL, which=NULL) {

  cluster <- NULL
  
  if (!is.null(h)) {
    cluster <- cutree(hcluster, h=h)
  } else if (!is.null(k)) {
    cluster <- cutree(hcluster, k=k)
  } 

  classlabel   <- as.factor(cluster)
    clusterreps  <- tapply(mols, classlabel, metfRag:::getMCSS)
    
    if (missing(which)) {
      which=seq(along=clusterreps)
    }

  m <- unique(cluster[hcluster$order])
  
  clusterreps[m[which]]
}

#' Chemical Clustering
#' 
#' Calculate the Maximum Common Substructure
#' Modified from stats::rect.hclust()
#' 
#' @aliases myimages.hclust
#' @param tree an object of the type produced by \code{hclust},
#' that was just plot()ed, and where the images should be overlaid
#' @param k,h Scalar. Cut the dendrogram such that either exactly
#' \code{k} clusters are produced or by cutting at height \code{h}.
#' (either k or h needs to be specified)
#' @param which,x A vector selecting the clusters around which a
#' rectangle should be drawn. \code{which} selects clusters by number
#' (from left to right in the tree), \code{x} selects clusters
#' containing the respective horizontal coordinates. Default is
#' \code{which = 1:k} (either x or which needs to be specified)
#' @param cluster Optional vector with cluster memberships as returned by
#' \code{cutree(hclust.obj, k = k)}, can be specified for efficiency if
#' already computed.
#' @param border Vector with border colors for the rectangles, NULL for none. Recycled if neccessary.
#' @param mols The a list of rCDK \code{mols}
#' 
#' @export


myimages.hclust <- function (tree, k = NULL, which = NULL, x = NULL, h = NULL,
                             cluster = NULL, mols=NULL, border = NULL) 
{
  if (length(h) > 1L | length(k) > 1L) 
    stop("'k' and 'h' must be a scalar")

  if (!is.null(h)) {
    if (!is.null(k)) 
      stop("specify exactly one of 'k' and 'h'")
    k <- min(which(rev(tree$height) < h))
    k <- max(k, 2)
  } else if (is.null(k)) {
    stop("specify exactly one of 'k' and 'h'")
  }
  
  if (k < 2 | k > length(tree$height)) {
    stop(gettextf("k must be between 2 and %d", length(tree$height)), 
         domain = NA)
  }
  
  if (is.null(cluster)) {
    cluster <- cutree(tree, k = k)
  }
  
  ## cutree returns classes sorted by data, we need classes
  ## as occurring in the tree (from left to right)  
  clustab <- table(cluster)[unique(cluster[tree$order])]
  m <- c(0, cumsum(clustab))
  if (!is.null(x)) {
    if (!is.null(which)) 
      stop("specify exactly one of 'which' and 'x'")
    which <- x
    for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
  } else if (is.null(which)) {
    which <- 1L:k
  }
  
  if (any(which > k)) 
    stop(gettextf("all elements of 'which' must be between 1 and %d", 
                  k), domain = NA)
  
  retval <- list()
  if (!is.null(border)) {
    border <- rep(border, length_out = length(which))
    for(n in seq_along(which)) {
      rect(m[which[n]]+0.66, par("usr")[3L],
           m[which[n]+1]+0.33, mean(rev(tree$height)[(k-1):k]),
           border = border[n])
      retval[[n]] <- which(cluster==as.integer(names(clustab)[which[n]]))
    }
  }
  
  retval <- list()
  for (n in seq_along(which)) {
    consensus <- getMCSS(mols[which(cluster == as.integer(names(clustab)[which[n]]))])
    xwidth <- par("usr")[2]-par("usr")[1]
    ywidth <- par("usr")[4]-par("usr")[3]
    aspect <- xwidth/ywidth
    
    xleft <- m[which[n]] + 0.66
    width <- m[which[n] + 1] + 0.33 - xleft
    ybottom <- par("usr")[3L]
    
    rasterImage(view.image.2d(consensus),
                xleft, 
                ybottom,
                xleft + width,        
                (ybottom + width)/aspect)
    
    retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
  }
  invisible(retval)
}

#' Chemical Clustering
#' 
#' Calculate the Maximum Common Substructure
#' Modified from stats::rect.hclust()
#' 
#' @aliases myimages.clustNumbers
#' @param tree an object of the type produced by \code{hclust},
#' that was just plot()ed, and where the images should be overlaid
#' @param k,h Scalar. Cut the dendrogram such that either exactly
#' \code{k} clusters are produced or by cutting at height \code{h}.
#' (either k or h needs to be specified)
#' @param which,x A vector selecting the clusters around which a
#' rectangle should be drawn. \code{which} selects clusters by number
#' (from left to right in the tree), \code{x} selects clusters
#' containing the respective horizontal coordinates. Default is
#' \code{which = 1:k} (either x or which needs to be specified)
#' @param cluster Optional vector with cluster memberships as returned by
#' \code{cutree(hclust.obj, k = k)}, can be specified for efficiency if
#' already computed.
#' @param border Vector with border colors for the rectangles, NULL for none. Recycled if neccessary.
#' 
#' @export


myimages.clustNumbers <- function (tree, k = NULL, which = NULL, x = NULL, h = NULL,
                             cluster = NULL, border = NULL) 
{
  if (length(h) > 1L | length(k) > 1L) 
    stop("'k' and 'h' must be a scalar")
  
  if (!is.null(h)) {
    if (!is.null(k)) 
      stop("specify exactly one of 'k' and 'h'")
    k <- min(length(tree$height)+1, which(rev(tree$height) < h)) # Can 
    k <- max(k, 2)
  } else if (is.null(k)) {
    stop("specify exactly one of 'k' and 'h'")
  }

  #if (k < 2 | k > length(tree$height)) {
  #  stop("k must be between 2 and %d", length(tree$height))
  #}
  
  if (is.null(cluster)) {
    cluster <- cutree(tree, k = k)
  }
  
  ## cutree returns classes sorted by data, we need classes
  ## as occurring in the tree (from left to right)  
  clustab <- table(cluster)[unique(cluster[tree$order])]
  m <- c(0, cumsum(clustab))
  if (!is.null(x)) {
    if (!is.null(which)) 
      stop("specify exactly one of 'which' and 'x'")
    which <- x
    for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
  } else if (is.null(which)) {
    which <- 1L:k
  }
  
  if (any(which > k)) 
    stop("all elements of 'which' must be between 1 and %d", k)
  
  retval <- list()
  if (!is.null(border)) {
    border <- rep(border, length_out = length(which))
    for(n in seq_along(which)) {    
      rect(m[which[n]]+0.66, par("usr")[3L],
           m[which[n]+1]+0.33, 
           ifelse(is.na(mean(rev(tree$height)[(k-1):k])), 0.5, mean(rev(tree$height)[(k-1):k])) ,
           border = border[n])

      cexonewidth <- strwidth(which[n])
      clusterwidth <- (m[which[n]+1]+0.33)-(m[which[n]]+0.66)
      
      text(x=(m[which[n]]+0.66)+clusterwidth/2,
           y=(ifelse(is.na(mean(rev(tree$height)[(k-1):k])), 0.5, mean(rev(tree$height)[(k-1):k]))-par("usr")[3L])/2,
        labels=which[n], col=border[n], cex=clusterwidth/cexonewidth)
      
      retval[[n]] <- which(cluster==as.integer(names(clustab)[which[n]]))
    }
  }
  
  invisible(retval)
}
#plot(cluster, hang=-1)
#myimages.clustNumbers(cluster, k=8, which=1:8, border=2)


myidentify <- function (x, FUN = NULL, N = 20, MAXCLUSTER = 20, DEV.FUN = NULL, ...) 
{
  cluster <- cutree(x, k = 2:MAXCLUSTER)
  retval <- list()
  oldk <- NULL
  oldx <- NULL
  DEV.x <- grDevices::dev.cur()
  for (n in 1L:N) {
    grDevices::dev.set(DEV.x)
    X <- locator(1)
    if (is.null(X)) 
      break
    k <- min(which(rev(x$height) < X$y), MAXCLUSTER)
    k <- max(k, 2)
    retval[[n]] <- unlist(myimages.hclust(x, k = k, x = X$x, 
                                          cluster = cluster[, k - 1]))
    if (!is.null(FUN)) {
      if (!is.null(DEV.FUN)) {
        grDevices::dev.set(DEV.FUN)
      }
      retval[[n]] <- FUN(retval[[n]], ...)
    }
    oldx <- X$x
    oldk <- k
  }
  grDevices::dev.set(DEV.x)
  invisible(retval)
}


# plot(cluster)
# myidentify(cluster)
