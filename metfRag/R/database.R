Sys.setlocale("LC_NUMERIC",'C');

db.search <- function(file, db, type, value)
{
  database <- c('kegg', 'pubchem');
  
  if (any(db == database) == TRUE)
  {
    func.getId    <- paste("db",db,"getId",sep=".");
    func.chemFile <- paste("db",db,"getMoleculeContainer",sep=".");
    
    db.ids <- do.call(func.getId, list(type, value));
    do.call(func.chemFile, list(file, db.ids));
  }
}

#KEGG
db.kegg.getId <- function(type, value)
{
  kegg.types  <- c('exact_mass', 'mol_weight', 'formula');
  kegg.loc    <- "http://rest.kegg.jp/find";
  kegg.db     <- "compound";  
  
  if (any(kegg.types == type) == TRUE)
  {
    if (kegg.types[1] == type)
    { value <- as.double(value); }
    
    kegg.url <- paste(kegg.loc,kegg.db,value,type,sep="/");
    kegg.data <- getURI(kegg.url);
    
    return(db.kegg.convertString(kegg.data));
  }
  
  return(NULL);
}

db.kegg.convertString <- function(idstring)
{
  ids   <- strsplit(idstring, "\n");
  size  <- length(ids[[1]]);
  kegg.id   <- NULL;
  tmp       <- NULL;
  converted <- sapply(ids[[1]],strsplit,"\t")
  
  for (i in 1:(size))
  {
    tmp <- strsplit(converted[[i]][1],":");
    kegg.id <- rbind(kegg.id, tmp[[1]][2]);
  }
  
  colnames(kegg.id) <- 'CPD';
  return(kegg.id);
}
  
db.kegg.getMoleculeContainer <- function(file=NULL, ids)
{
  if (length(ids) == 0)
  { return(FALSE); }  
  
  ids <- paste("cpd:",ids,sep="");
  
  kegg.loc  <- "http://rest.kegg.jp";
  kegg.op   <- 'get';
  kegg.filetype <- 'mol';
  kegg.id <- paste(ids, collapse="+");        
  
  kegg.url     <- paste(kegg.loc,kegg.op,kegg.id,kegg.filetype,sep="/"); 
  
  if (missing(file) == TRUE)
  {
    if (url.exists(kegg.url) == TRUE)
    { return(load.molecules(kegg.url)); }    
  }
  else
  {
    download.file(url, file);
    
    if (file.exists(file))
    { return(TRUE); }    
  }
  
  return(FALSE);
}

#PubChem
db.pubchem.getId <- function(type=NULL, value=NULL)
{
  pc.loc <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";
  pc.query <- "db=pccompound&term=";
  pc.type  <- c('formula','mimass','exmass','cid');
  
  if (any(pc.type == type) == TRUE && missing(value) != TRUE)
  {
    if (type == pc.type[2] || type == pc.type[3])
    {
      value <- as.double(value);
      pc.mass <- paste(100,
                       "[",toupper(type),"]:",
                       value,
                       "[",toupper(type),"]",
                       sep="");
    }
    else if (type == pc.type[4])
    {
      pc.mass <- paste("[",toupper(type),"]",
                       value,
                       "[",toupper(type),"]",
                       sep="");      
    }
    else
    { pc.mass <- value; }
    
    pc.url <- paste(pc.loc, pc.query, pc.mass, sep="");
    
    pc.data <- getURI(pc.url);
    return(db.pubchem.XMLToList(pc.data));
  }
  
  return(NULL);
}  

db.pubchem.XMLToList <- function(xml)
{
  pc.data <- xmlToList(xmlParse(xml));
  
  if (is.null(pc.data$WarningList) == TRUE)
  { 
    pc.data <- as.matrix(pc.data$IdList);
    colnames(pc.data) <- 'CID';
    rownames(pc.data) <- NULL;
    
    return(pc.data);  
  }
  
  return(NULL);
}

db.pubchem.getMoleculeContainer <- function(file=NULL, ids)
{
  pc.loc      <- "http://pubchem.ncbi.nlm.nih.gov/rest/pug";
  pc.db       <- "compound"
  pc.type     <- "cid";
  pc.filetype <- "SDF";
  pc.ids <- paste(ids,collapse=",");
  
  pc.url <- paste(pc.loc,pc.db,pc.type,pc.ids,pc.filetype,sep="/");
  
  if (missing(file) == TRUE)
  {
    if (url.exists(pc.url) == TRUE)
    {
      return(load.molecules(pc.url));
    }    
  }
  else
  {
    download.file(pc.url, file);
    
    if (file.exists(file))
    { return(TRUE); }    
  }  
  
  return(FALSE);
}