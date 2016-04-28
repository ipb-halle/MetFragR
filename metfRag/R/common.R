common.lib.findMolecule <- function(mol, property, seek.value, split=NULL)
{
  mol.value <- get.property(mol, property);
  #stuff <- names(get.properties(mol));
  #stop(stuff)
  
  if (is.null(split) == TRUE && mol.value == seek.value)
  {
    return(TRUE);
  }
  else if (missing(split) == FALSE && 
           is.null(split$sep) == FALSE && 
           is.null(split$pos) == FALSE)
  {
      mol.value <- strsplit(mol.value,split=split$sep)[[1]][split$pos];
      
      if (mol.value == seek.value)
      { return(TRUE); }
  }
  
  return(FALSE);
}

common.lib.lookup <- function(needle, haystack, link, split=NULL)
{
  if (missing(needle) == TRUE || 
      missing(haystack) == TRUE || 
      length(haystack) == 0)
  { return(FALSE); }
  
  occur <- c();
  
  if (is.null(split) == FALSE && 
      is.null(split$sep) == FALSE && 
      is.null(split$pos) == FALSE)
  {
    mol.value.needle <- strsplit(get.property(needle, link),
                                 split=split$sep)[[1]][split$pos];    
    
    for (i in (1:length(haystack)))
    {
      mol.value.hay <- strsplit(get.property(haystack[[i]], link),
                                split=split$sep)[[1]][split$pos];      
      
      if (mol.value.needle == mol.value.hay)
      { occur <- c(occur, i); }
    }
  }
  else
  {
    for (i in (1:length(haystack)))
    {
      if (get.property(needle, link) == get.property(haystack[[i]], link))
      { occur <- c(occur, i); }
    }
  }
  
  return (occur);
}

common.lib.lookupFirst <- function(needle, haystack, link, split=NULL)
{
  if (missing(needle) == TRUE || 
      missing(haystack) == TRUE || 
      length(haystack) == 0)
  { return(FALSE); }

  if (is.null(split) == FALSE && 
        is.null(split$sep) == FALSE && 
        is.null(split$pos) == FALSE)
  {
    mol.value.needle <- strsplit(get.property(needle, link),
                                 split=split$sep)[[1]][split$pos];    
    
    for (i in (1:length(haystack)))
    {
      mol.value.hay <- strsplit(get.property(haystack[[i]], link),
                                split=split$sep)[[1]][split$pos];      
      
      if (mol.value.needle == mol.value.hay)
      { return(TRUE); }
    }
  }
  else
  {
    for (i in (1:length(haystack)))
    {
      if ((get.property(needle, link) == get.property(haystack[[i]], link)))
      { return(TRUE); }
    }
  }
  
  return(FALSE);
}

lib.CommonSet <- function(set)
{
  set.result <- c();
  
  for (i in 1:length(set))
  { set.result <<- c(set.result, names(set[[i]])); }

  set.result <- intersect(set.result, set.result);

  return(set.result);
}

comm.lib.showLinkOptions <- function(set.a, set.b)
{
  if (missing(set.a) == TRUE || missing(set.b) == TRUE)
  { return(FALSE); }
  
  set.a.prop <- sapply(set.a, get.properties);
  set.b.prop <- sapply(set.b, get.properties);
  
  set.inter <- intersect(
    lib.CommonSet(set.a.prop),
    lib.CommonSet(set.b.prop)
  );
    
  return(set.inter);
}

comm.lib.showNumberOptions <- function(set)
{
  if (missing(set) == TRUE)
  { return(FALSE); }
  
  keys <- comm.lib.showLinkOptions(set,set);
  name <- c();
  
  for (i in (1:length(keys)))
  {
    if (is.numeric(get.property(set[[1]],keys[i])))
    { name <- c(name, keys[i])}
    else
    {
      x <- as.numeric(get.property(set[[1]],keys[i]));
      
      if (is.na(x) == FALSE)
      {name <- c(name, keys[i]);}
    }
  }
  
  return(name);
}
