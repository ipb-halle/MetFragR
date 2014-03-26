Sys.setlocale("LC_NUMERIC",'C');

db.search <- function(file, db, type, value)
{
  database <- c('kegg', 'pubchem');
  
  if (any(db == database) == TRUE)
  {
    func.getId    <- paste("db",db,"getId",sep=".");
    func.chemFile <- paste("db",db,"MoleculeToFile",sep=".");
    
    db.ids <- do.call(func.getId, list(type, value));
    do.call(func.chemFile, list(file, db.ids));
  }
}

#KEGG
db.kegg.getId <- function(type, value)
{
  kegg.types  <- c('exact_mass', 'formula');
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
  ids   <- str_split(idstring, "\n");
  size  <- length(ids[[1]]);
  kegg.id   <- NULL;  
  converted <- sapply(ids[[1]],str_split,"\t")
  
  for (i in 1:(size-1))
  {
    kegg.id <- rbind(kegg.id, converted[[i]][1]);
  }
  
  colnames(kegg.id) <- 'CPD';
  return(kegg.id);
}
  
db.kegg.MoleculeToFile <- function(file, ids)
{
  if (length(ids) == 0)
  { return(FALSE); }  
  
  kegg.loc  <- "http://rest.kegg.jp";
  kegg.op   <- 'get';
  kegg.filetype <- 'mol';
  kegg.id <- paste(ids, collapse="+");        
  
  url     <- paste(kegg.loc,kegg.op,kegg.id,kegg.filetype,sep="/"); 
  download.file(url, file);
  
  if (file.exists(file))
  { return(TRUE); }
  
  return(FALSE);
}

#PubChem
db.pubchem.getId <- function(type, value=NULL)
{
  pc.loc <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";
  pc.query <- "db=pccompound&term=";
  pc.type  <- c('formula', 'mimass','exmass');
  
  if (any(pc.type == type) == TRUE && missing(value) != TRUE)
  {
    if (type == pc.type[2] || type == pc.type[3])
    { value = as.double(value); }
    
    if (type != pc.type[1])
    {
      pc.mass <- paste(100,
                       "[",toupper(type),"]:",
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
  pc.data <- as.matrix(pc.data$IdList);
  colnames(pc.data) <- 'CID';
  rownames(pc.data) <- NULL;
  
  return(pc.data);
}

db.pubchem.MoleculeToFile <- function(file, ids)
{
  pc.loc      <- "http://pubchem.ncbi.nlm.nih.gov/rest/pug";
  pc.db       <- "compound"
  pc.type     <- "cid";
  pc.filetype <- "SDF";
  pc.ids <- paste(ids,collapse=",");
  
  pc.url <- paste(pc.loc,pc.db,pc.type,pc.ids,pc.filetype,sep="/");
  download.file(pc.url, file);
  
  if (file.exists(file))
  { return(TRUE); }  
  
  return(FALSE);
}