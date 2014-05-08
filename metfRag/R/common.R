common.lib.lookup <- function(needle, haystack, link)
{
  if (missing(needle) == TRUE || 
      missing(haystack) == TRUE || 
      length(haystack) == 0)
  { return(FALSE); }
  
  occur <- c();
  
  for (i in (1:length(haystack)))
  {
    if ((get.property(needle, link) == get.property(haystack[[i]], link)))
    { 
      occur <- c(occur, i); 
    }
  }
  
  return (occur);
}

common.lib.lookupFirst <- function(needle, haystack, link)
{
  if (missing(needle) == TRUE || 
      missing(haystack) == TRUE || 
      length(haystack) == 0)
  { return(FALSE); }
  
  for (i in (1:length(haystack)))
  {
    if ((get.property(needle, link) == get.property(haystack[[i]], link)))
    { 
      return(TRUE);
    }
  }
  
  return(FALSE);
}

comm.lib.showLinkOptions <- function(set)
{
  if (missing(set) == TRUE)
  { return(FALSE); }
  
  name <- names(get.properties(set[[1]]));
  return(name);
}