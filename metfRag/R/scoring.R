scoring.getMolValues <- function(sorting, mols)
{  
  sorter <- c();

  for (i in (1:length(mols)))
  { sorter <- c(sorter, as.numeric(get.property(mols[[i]], sorting))); }
  
  return(sorter);
}

scoring.applyFunc <- function(sorting, mols)
{
  column <- sapply(sorting[c(1:(length(sorting)-1))], scoring.getMolValues, mols);
  molrow <- c();
  
  for (i in (1:nrow(column)))
  {
    molrow <- c(molrow, do.call(
      sorting[[length(sorting)]], 
      as.list(column[i, 1:ncol(column)])
    ));
  }
  
  return(molrow);
}

scoring.sortRangOrder <- function(molMatrix)
{
  common.columns <- as.numeric(colnames(molMatrix));
  mol.columns <- NULL;
  
  for(i in (1:(length(common.columns))))
  { mol.columns <- cbind(mol.columns, molMatrix[,common.columns[i]]); }
  
  return(mol.columns);
}

is.integer0 <- function(val)
{
  if (is.integer(val) == TRUE && length(val) == 0L)
  {
    return(TRUE);
  }
  
  return(FALSE);
}

scoring.orderContainer <- function(sorting, mols)
{
  molprop <- c();
  
  sort.types  <- sapply(sorting, class);
  pos.list    <- which(sort.types=='list');
  pos.char    <- which(sort.types=='character');
  col.list    <- NULL;
  col.char    <- NULL;
  
  if (is.integer0(pos.char) == FALSE)
  {
    col.char <- sapply(sorting[pos.char], scoring.getMolValues, mols);
    colnames(col.char) <- pos.char;    
  }
  
  if (is.integer0(pos.list) == FALSE)
  {
    col.list <- sapply(sorting[pos.list], scoring.applyFunc, mols);
    colnames(col.list) <- pos.list;
  }
  
  mol.list <- cbind(col.char, col.list);
  mol.list <- as.data.frame(scoring.sortRangOrder(mol.list));
  
  mol.list <- na.omit(
    mol.list[order(mol.list[,1:length(mol.list)], decreasing=TRUE),]);
  
  return(mol.list);
}

scoring.calcCategories <- function(target.position, rank.frame)
{
  rank.frame <- as.matrix(rank.frame);
  colnames(rank.frame) <- NULL;
  rownames(rank.frame) <- NULL;
  
  rank <- list(BC=0, WC=0, EC=0, TC=0); 
  
  for (i in (1:nrow(rank.frame)))
  {
    if (i == target.position)
    { next; }
    
    better <- pmax(rank.frame[target.position,], rank.frame[i,]);
    worst  <- pmin(rank.frame[target.position,], rank.frame[i,]);
    
    if (isTRUE(all.equal(rank.frame[target.position,], rank.frame[i,])) == TRUE)
    { rank$EC = rank$EC + 1; }
    else if(isTRUE(all.equal(rank.frame[i,],worst)) == TRUE)
    { rank$WC = rank$WC + 1; }
    else if(isTRUE(all.equal(rank.frame[i,],better)) == TRUE) 
    { rank$BC = rank$BC + 1; }
  }
  
  rank$TC = nrow(rank.frame);

  return(rank);
}

scoring.calcMolParameter <- function(row.pos, mol.pos, orderedContainer)
{
  rank <- scoring.calcCategories(row.pos, orderedContainer);
  rank <- c(
    rank,
    pessimistic = rank$BC + rank$EC + 1,
    optimistic = rank$BC + 1,
    rrp = 0.5 * (1- ((rank$BC-rank$WC)/(rank$TC - 1))),
    molecule_position = mol.pos
  )
}

scoring.getRanks <- function(mols, sorting, condition, split=NULL)
{
  if (missing(mols) == TRUE)
  { return(FALSE); }  
    
  if (missing(sorting) == TRUE)
  {
    	warning("Choose valid sorting parameter.")
	return(FALSE);
  }
  
  if (missing(condition) == TRUE)
  {
	 warning("Choose valid condition parameter.")
    	return(FALSE);
  }
  
  orderedContainer  <- scoring.orderContainer(sorting, mols);
  if(is.null(dim(orderedContainer))) {orderedContainer=data.frame(V1=orderedContainer)}
  mol.pos <- sapply(mols, common.lib.findMolecule, 
                    names(condition), condition[[1]], split);
  
  mol.pos           <- which(mol.pos);
  pos.rows          <- match(mol.pos, as.numeric(rownames(orderedContainer)));
  
  rank <- mapply(scoring.calcMolParameter, pos.rows, mol.pos, 
                 MoreArgs=list(orderedContainer));
  return(rank);
}
