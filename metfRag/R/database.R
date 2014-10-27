Sys.setlocale("LC_NUMERIC",'C');

#Library
db.lib.matchFunc <- function(names.seek, names.func)
{
  is.matched <- match(names.func, names.seek);
  idx <- NULL;
  
  for (i in 1:length(is.matched))
  {
    if (is.na(is.matched[i]) == FALSE)
    {
      idx <- c(idx,i);
    }
  }
  
  if (is.null(idx) == FALSE)
  { return(idx); }
    
  return(NULL);
}

db.lib.calcMass <- function(seek)
{
  value <- as.double(seek[1]);
  
  if (is.null(seek[2]) == FALSE)
  {
    range <- as.double(seek[2]);
  }
  else
  { range <- as.double(0); }
   
  seek <- list(
    val = value,
    lbound=(value - range),
    ubound=(value + range)
  );
  
  return(seek);
}

#ChemSpider
db.chemspider.buildURI <- function(func, values)
{
  cs.loc <- "http://www.chemspider.com/MassSpecAPI.asmx";
  cs.func <- func;
  cs.main <- paste(cs.loc, cs.func, sep="/");
  cs.values <- paste(names(values), values, sep="=");
  
  if(length(values) > 1)
  { cs.values <- paste(cs.values,collapse="&"); }
  
  cs.url <- paste(cs.main, cs.values, sep="?");
  
  return(cs.url);
}

db.chemspider.XMLToList <- function(xml)
{
  xmlList <- xmlToList(xml);
  xmlList <- as.matrix(xmlList);
  colnames(xmlList) <- "CSID";
  rownames(xmlList) <- NULL;
  
  return(xmlList);
}

db.chemspider.getId <- function(seek=NULL)
{
  if (is.null(seek) == TRUE)
  { return(FALSE); }
  
  cs.seek.names <- names(seek);
  cs.seek.funcs <- list(mass = "SearchByMass2", 
                        formula = "SearchByFormula2");
  
  idx <- db.lib.matchFunc(cs.seek.names, names(cs.seek.funcs));
  
  if (is.null(seek$range) == TRUE && is.null(seek$mass) == FALSE)
  { seek = c(seek, range=as.double(0)); }
  
  if (is.null(idx) == FALSE)
  { cs.url <- db.chemspider.buildURI(cs.seek.funcs[idx], seek); }
  
  cs.data <- getURI(cs.url);
  cs.data <- db.chemspider.XMLToList(cs.data);
  
  return(cs.data);
}

db.chemspider.getMoleculeContainer <- function(ids, token, cal3d=FALSE)
{
  cs.func <- list(mol = "GetRecordMol");
  cs.idList <- NULL;
  
  for (i in 1:length(ids))
  {
    cs.idList <- c(cs.idList,
                   db.chemspider.buildURI(cs.func, 
                                          list(csid=ids[[i]][1], 
                                               calc3d=cal3d, 
                                               token=token)));
  } 
  return(cs.idList)
  xmlLists <- getURIAsynchronous(cs.idList);
  x <- db.chemspider.XMLToList(xmlLists[1])
  #x <- load.molecules(x)
  #x <- lapply(xmlLists, db.chemspider.XMLToList);
  
  return(x);
}

#KEGG
db.kegg.getQuery <- function(seek, types)
{
  names.seek <- names(seek);
  idx <- db.lib.matchFunc(names.seek, names(types));
  
  if (is.null(idx) == FALSE)
  {    
    seek <- list(str=seek[[1]], idx=idx);
    
    if (is.null(seek$mass) == FALSE)
    { 
      seek <- db.lib.calcMass(seek);
      seek <- list(str=paste(seek$lbound,seek$ubound,sep="-"), idx=idx);
    }  
    else if (is.null(seek$mw) == FALSE)
    {
      seek <- db.lib.calcMass(seek);
      seek <- list(str=paste(seek$lbound,seek$ubound,sep="-"), idx=idx);      
    }
    
    return(seek);
  }
  
  return(FALSE);
}

db.kegg.getId <- function(seek=NULL)
{
  if(is.null(seek) == TRUE)
  { return(FALSE); }
  
  kegg.types <- list(mass = 'exact_mass', 
                     mw = 'mol_weight', 
                     formula = 'formula');
  
  seek <- db.kegg.getQuery(seek, kegg.types);

  kegg.loc   <- "http://rest.kegg.jp/find";
  kegg.db    <- "compound";  
    
  kegg.url <- paste(kegg.loc, kegg.db, seek$str, kegg.types[seek$idx], sep="/");
  
  if (url.exists(kegg.url) == TRUE)
  { 
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

db.kegg.annotateMolecule <- function(container, id)
{  
  file <- keggGet(dbentries=c(id));
  skip.properties <- c("BOND", "ATOM", "BRITE", "ENZYME" ,"DBLINKS", 
                       "COMMENT");
  
  for (i in ( 1:length(names(file[[1]])) ))
  {
    mol.name <- names(file[[1]][i]);
    mol.value <- file[[1]][[i]];
    
    if (mol.name == "ENTRY")
    { mol.value <- strsplit(mol.value, split="\n")$Compound; }
    else if (any(mol.name == skip.properties) == TRUE)
    { next; }
    else if(any(mol.name == c("REACTION","NAME")) == TRUE)
    { mol.value <- paste(mol.value, collapse=" "); }
    else if(mol.name == "PATHWAY")
    {
      for (j in (1:length(file[[1]]$PATHWAY)))
      {
        mol.name <- paste("KEGG", "PATHWAY", 
                          names(file[[1]]$PATHWAY[j]), sep="_");
        mol.val <- file[[1]]$PATHWAY[[j]];
        set.property(container, mol.name, mol.val);
      }
      
      next;
    }
      
    mol.name <- paste("KEGG", names(file[[1]][i]), sep="_");
    set.property(container, mol.name, mol.value);
  }
  
  return(container);
}

db.kegg.getMoleculeContainer<- function(ids=NULL)
{
  if (length(ids) == 0 || is.null(ids) == TRUE)
  { return(FALSE); }  
  
  ids           <- paste("cpd:",ids,sep="");
  kegg.loc      <- "http://rest.kegg.jp";
  kegg.op       <- 'get';
  kegg.filetype <- 'mol';
  kegg.url      <- NULL;
  
  kegg.url <- sapply(
    ids,
    function(id, type, loc, op) 
      paste(loc, op, id, type, sep="/"),
    kegg.filetype,
    kegg.loc,
    kegg.op
  );
  
  molecules <- sapply(kegg.url, load.molecules);
  molecules <- mapply(db.kegg.annotateMolecule, molecules, names(molecules));
  names(molecules) <- NULL;
  molecules <- lapply(molecules, container.identicalProperties);
  
  return(molecules);
}

#PubChem
db.pubchem.getId <- function(seek=NULL)
{  
  pc.loc <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";
  pc.query <- "db=pccompound&term=";
  pc.maxresult <- paste("retmax", 100, sep="=");
  
  pc.type  <- list(formula='formula',
                   mass='exmass',
                   cid='cid');
  
  pc.idx <- db.lib.matchFunc(names(seek), names(pc.type));
  
  if (is.null(seek) == TRUE || length(seek) > 2 || is.null(pc.idx) == TRUE )
  { return(FALSE); }
  
  if (is.null(seek$mass) == FALSE)
  {
    seek <- db.lib.calcMass(seek);

    pc.mass <- paste(seek$lbound,
                     "[",toupper(pc.type[pc.idx]),"]:",
                     seek$ubound,
                     "[",toupper(pc.type[pc.idx]),"]",
                     sep="");
  }
  else if (is.null(seek$cid) == FALSE)
  {
    pc.mass <- paste("[",toupper(pc.type[pc.idx]),"]",
                     seek$cid,
                     "[",toupper(pc.type[pc.idx]),"]",
                     sep="");      
  }
  else
  { pc.mass <- seek[1]; }
    
  pc.url <- paste(pc.loc, pc.query, pc.mass, sep="");
  pc.url <- paste(pc.url, pc.maxresult, sep="&")
  pc.data <- getURI(pc.url);
  
  return(db.pubchem.XMLToList(pc.data));
}  

db.pubchem.XMLToList <- function(xml)
{
  pc.data <- xmlToList(xmlParse(xml));
  
  if (is.null(pc.data$WarningList) == TRUE && pc.data$Count != '0')
  { 
    pc.data <- as.matrix(pc.data$IdList);
    colnames(pc.data) <- 'CID';
    rownames(pc.data) <- NULL;
    
    return(pc.data);  
  }
  
  warning(pc.data$WarningList);
  return(NULL);
}

db.pubchem.getMoleculeContainer <- function(ids=NULL)
{
  if (is.null(ids) == TRUE)
  { return(FALSE); }
  
  pc.loc      <- "http://pubchem.ncbi.nlm.nih.gov/rest/pug";
  pc.db       <- "compound"
  pc.type     <- "cid";
  pc.filetype <- "SDF";
  pc.ids <- paste(ids,collapse=",");
  
  pc.url <- paste(pc.loc,pc.db,pc.type,pc.ids,pc.filetype,sep="/");
  
  if (url.exists(pc.url) == TRUE)
  { 
    molecules <- load.molecules(pc.url);
    molecules <- lapply(molecules, container.identicalProperties);
    return(molecules); 
  }
  
  return(FALSE);
}