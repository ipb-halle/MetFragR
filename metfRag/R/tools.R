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
  return (sapply(mols, function(x) as.numeric(get.property(x, "NoPeaksExplained"))))
}

#' MetFrag result tools
#' 
#' Small tools to extract information from MetFrag or MetFusion results
#' 
#' @aliases getDatabaseIDs
#' @param mols A list of rCDK \code{molecules} with MetFrag or MetFusion results
#' @param idprop The name of the property of the molecules where the database ID is kept
#' @export
getDatabaseIDs <- function(mols, idprop="DatabaseID") {
  return (sapply(1:length(mols), function(i) {
    id <- get.property(mols[[i]], idprop)
    if (is.na(id) || id == 'NA')
      return (i)
    else
      return (get.property(mols[[i]], idprop))
  }))
}
