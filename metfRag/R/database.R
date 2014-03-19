Sys.setlocale("LC_NUMERIC",'C');

db.searchToFile <- function(file, db, type, value)
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
  kegg.types <- c('exact_mass', 'entry', 'name', 'formula');
  
  if (any(kegg.types == type) == TRUE)
  {
    if (kegg.types[1] == type)
    { value <- as.double(value); }
    
    kegg.id <- names(keggFind("compound", value, type));
    return(kegg.id);
  }
  
  return(c());
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
db.pubchem.getId <- function(type, value)
{
  pc.loc <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";
  pc.query <- "db=pccompound&term=";
  pc.type  <- c('formula', 'mimass','exmass', 'cid');
  
  if (any(pc.type == type) == TRUE)
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
    { pc.mass <- NULL; }
    
    pc.url <- paste(pc.loc, pc.query, pc.mass, sep="");
    pc.data <- getURI(pc.url);
    
    return(db.pubchem.XMLToList(pc.data))
  }
  
  return(c());
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