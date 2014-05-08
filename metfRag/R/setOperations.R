container.union <- function(set.a, set.b, link)
{  
  if (missing(link) == TRUE)
  {
    print(comm.lib.showLinkOptions(set.a));
    stop("Please choose a link from above. You could use the
         function: comm.lib.showLinkOptions(set)");
  }
  
  sets <- c(set.a, set.b);
  double <- lapply(sets, comm.lib.lookup, sets, link);
  result <- c();
  
  for (i in (1:length(double)))
  {
    if (length(double[[i]]) == 1)
    {
      result <- c(result, sets[double[[i]][i]]);
    }
    else if (length(double[[i]]) > 1)
    {
      if (i == double[[i]][1])
      { result <- c(result, sets[double[[i]][i]]); }
    }
  }
  
  return(result);
}

container.intersect <- function(set.a, set.b, link)
{  
  if (missing(link) == TRUE)
  {
    print(comm.lib.showLinkOptions(set.a));
    stop("Please choose a link from above. You could use the
         function: comm.lib.showLinkOptions(set)");
  }
  
  double <- lapply(set.a, comm.lib.lookup, set.b, link);
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
    print(comm.lib.showLinkOptions(set.a));
    stop("Please choose a link from above. You could use the
         function: comm.lib.showLinkOptions(set)");
  }
  
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
    print(comm.lib.showLinkOptions(set.a));
    stop("Please choose a link from above. You could use the
         function: comm.lib.showLinkOptions(set)");
  }
  
  result <- c();
  result <- c(result, container.asymmetric.difference(set.a, set.b, link));
  result <- c(result, container.asymmetric.difference(set.b, set.a, link));
  
  return(result);
}