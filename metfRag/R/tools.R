#' MetFrag result tools
#' 
#' Small tools to extract information from MetFrag or MetFusion results
#' 
#' @aliases getScores
#' @param mols A list of rCDK \code{molecules} with MetFrag or MetFusion results
#' @param scoreprop The name of the property of the molecules where the score is kept
#' @export

getScores <- function(mols, scoreprop="Score"){
  return (sapply(mols, function(x) as.numeric(get.property(x, scoreprop))))
}

#' MetFrag result tools
#' 
#' Small tools to extract information from MetFrag or MetFusion results
#' 
#' @aliases getPeaksExplained 
#' @param mols A list of rCDK \code{molecules} with MetFrag or MetFusion results
#' @export
getPeaksExplained <- function(mols){
  return (sapply(mols, function(x) as.numeric(get.property(x, "NoExplPeaks"))))
}

#' MetFrag result tools
#' 
#' Small tools to extract information from MetFrag or MetFusion results
#' 
#' @aliases getDatabaseIDs
#' @param mols A list of rCDK \code{molecules} with MetFrag or MetFusion results
#' @param idprop The name of the property of the molecules where the database ID is kept
#' @export
getDatabaseIDs <- function(mols, idprop="Identifier") {
  return (sapply(1:length(mols), function(i) {
    id <- get.property(mols[[i]], idprop)
    if (is.na(id) || id == 'NA')
      return (i)
    else
      return (get.property(mols[[i]], idprop))
  }))
}

#' Plotting a molecule
#' 
#' rCDK provides either functions to plot into a Java window,
#' or into a rasterImage object. plotMol
#' is a high-level plotting function to plot into an R device.
#' 
#' @aliases plotMol
#' @param mol An rCDK \code{mol} object
#' @param smiles A SMILES string, ignored if mol is supplied
#' @param width,height Number of pixels for molecule image in x,y
#' @param watermark A string written across the image as watermark
#' @export
plotMol <- function(mol=NULL, smiles=NULL, width = 200, height = 200, watermark=NULL) {
  if (missing(mol)) {
    mol <- parse.smiles(smiles)[[1]]
  }
  #  img <- view.image.2d(mol, width = 1024, height = 1024)
  img <- view.image.2d(mol, width = width, height = height)
  plot.new()
  plot.window(xlim=c(0,1), ylim=c(0,1))
  if (!is.null(watermark)) {
    text(x=0.5, y=0.5, label=watermark, 
         cex=0.55/strwidth(watermark), 
         col=rgb(0.75, 1,0.75))
  }
  rasterImage(alpha2image(img), 0,0, 1,1)
}

#' Image manipulation
#' 
#' This function adds an alpha channel to an image
#' 
#' @aliases alpha2image
#' @param img A rasterImage
#' @param threshold An intensity threshold, all pixels brighter than this 
#' will be set to transparent
#' @export
alpha2image <- function(img, threshold=0.5) {
  a2 <- as.array(img)
  th <- 0.5 ## Threshhold
  alpha <- a2[,,1]<th | a2[,,2]<th | a2[,,3]<th
  
  d <- dim(a2)
  d[3] <- 4
  a3 <- array(dim=d)
  
  a3[,,1:3] <- a2
  a3[,,4] <- alpha
  
  img <- as.raster(a3)  
}
  

#' Chemical Clustering
#' 
#' Calculate the Maximum Common Substructure for a set of molecules.
#' This is a star-alignment, iterating through all structures
#' and creating the consensus structure as Maximum Common Substructure 
#' between the consensus and the next molecule in the list.
#' 
#' @aliases getMCSS
#' @usage getMCSS(mols)
#' @param mols The a list of rCDK \code{mols}
#' @author Steffen Neumann (\email{sneumann@@ipb-halle.de})
#' @examples 
#'        library(rcdk)
#'        smiles <- c('CCC', 'CCN', 'CCN(C)(C)',
#'                    'c1ccccc1Cc1ccccc1',
#'                    'C1CCC1CC(CN(C)(C))CC(=O)CC')
#'        mols <- parse.smiles(smiles)
#'        dummy <- getMCSS(mols)
#' 
#' @export

getMCSS <- function(mols) {
  l <- length(mols)
  if (l<1)
    stop("No MCSS for less than one molecule!")
  
  if (l==1)
    return (mols[[1]])
  
  consensus <- mols[[1]]  
  for (i in 2:l) {
    consensus <- get.mcs(consensus, mols[[i]])
  }
  return (consensus)  
}
