container.identicalProperties <- function(mol)
{
  mol.prop <- get.properties(mol);
  mol.names <- names(mol.prop);
  
  for (i in (1:length(mol.names)))
  {
    if ( length(grep("^cdk:.*",mol.names[i])) == 1 )
    { next; }
    
    spl.pc <- strsplit(mol.names[i],split="PUBCHEM_")[[1]][2];
    spl.kegg <- strsplit(mol.names[i],split="KEGG_")[[1]][2];
    
    if (is.na(spl.pc) == FALSE)
    {     
      remove.property(mol, mol.names[i]);
      mol.names[i] <- spl.pc; 
      set.property(mol, mol.names[i], mol.prop[[i]]);
    }
    else if (is.na(spl.kegg) == FALSE)
    {     
      remove.property(mol, mol.names[i]);
      mol.names[i] <- spl.kegg; 
      set.property(mol, mol.names[i], mol.prop[[i]]);
    }
  }
  
  return(mol);
}

container.union <- function(set.a, set.b, link)
{    
  sets <- c(set.a, set.b);
  
  if (missing(link) == TRUE)
  {
    print(comm.lib.showLinkOptions(sets,sets));
    warning("Please choose a link from above. You could use the
         function: comm.lib.showLinkOptions(set.a,set.b)");
    return(FALSE);
  }
  
  link <- toupper(link);
  double <- lapply(sets, common.lib.lookup, sets, link);
  result <- c();
  
  for (i in (1:length(double)))
  {
    if (length(double[[i]]) == 1)
    {
      result <- c(result, sets[double[[i]][i]]);
    }
    else if (length(double[[i]]) > 1)
    {
      if (i == double[[i]][1] && is.na(double[[i]][i]) == FALSE)
      { 
        result <- c(result, sets[double[[i]][i]]); 
      }
    }
  }
  
  return(result);
}

container.intersect <- function(set.a, set.b, link)
{   
  if (missing(link) == TRUE)
  {
    print(comm.lib.showLinkOptions(set.a,set.b));
    warning("Please choose a link from above. You could use the
         function: comm.lib.showLinkOptions(set.a,set.b)");
    return(FALSE);
  }
  
  link <- toupper(link);
  double <- lapply(set.a, common.lib.lookup, set.b, link);
  result <- c();
  g <- 0;
  
  for (i in (1:length(double)))
  {
    if (length(double[[i]]) >= 1)
    {
      if (length(result) > 0 && 
          common.lib.lookupFirst(set.a[[i]], result, link) == FALSE)
      {
        result <- c(result, set.a[[i]]);
      }
      else if (length(result) == 0)
      {
        result <- c(result, set.a[[i]]);
      }
    }
    
  }
  
  return(result);
}

container.asymmetric.difference <- function(set.a, set.b, link)
{  
  if (missing(link) == TRUE)
  {
    print(comm.lib.showLinkOptions(set.a,set.b));
    warning("Please choose a link from above. You could use the
         function: comm.lib.showLinkOptions(set.a,set.b)");
    return(FALSE);
  }
  
  link <- toupper(link);
  double <- lapply(set.a, common.lib.lookupFirst, set.b, link);
  result <- c();
  
  for (i in (1:length(double)))
  {
    if (double[[i]] == FALSE)
    {
      if (length(result) == 0)
      {
        result <- c(result, set.a[[i]]);
      }
      else if (common.lib.lookupFirst(set.a[[i]], result, link) == FALSE)
      {
        result <- c(result, set.a[[i]]);        
      }
    }
  }
  
  return(result);  
}

container.symmetric.difference <- function(set.a, set.b, link)
{  
  if (missing(link) == TRUE)
  {
    print(comm.lib.showLinkOptions(set.a,set.b));
    warning("Please choose a link from above. You could use the
         function: comm.lib.showLinkOptions(set.a,set.b)");
    return(FALSE);
  }
  
  link <- toupper(link);
  result <- c();
  result <- c(result, container.asymmetric.difference(set.a, set.b, link));
  result <- c(result, container.asymmetric.difference(set.b, set.a, link));
  
  return(result);
}