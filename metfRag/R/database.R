library("rcdk");
library("SSOAP");
library("KEGGREST");
library("RCurl");

Sys.setlocale("LC_NUMERIC",'C');

#KEGG
db.kegg.getIdByMass <- function(mass)
{
  mass      <- as.double(mass);
  kegg.id   <- names(keggFind("compound",mass,"exact_mass"));
  
  return(kegg.id);
}

db.kegg.getIdByFormula <- function(formula)
{
  mass      <- as.double(mass);
  kegg.id   <- names(keggFind("compound",formula,"formula"));
  
  return(kegg.id);
}

db.kegg.getMolecule <- function(kegg.id)
{
  kegg.mol <- NULL;
  
  for (i in (1:length(kegg.id)))
  {
    kegg.url    <- paste("http://rest.kegg.jp/get/",kegg.id[i],"/mol",sep="");
    kegg.mol    <- c(kegg.mol, load.molecules(molfiles=kegg.url));    
  }
  
  return(kegg.mol);
}

#PubChem
db.pubchem.getIdByMass <- function(mass)
{  
  pc.loc    <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";
  pc.query  <- "db=pccompound&term=100";
  #[EXMASS]
  pc.mass   <- paste("[MIMASS]:", mass, "[MIMASS]", sep="");
  pc.url    <- paste(pc.loc,pc.query,pc.mass,sep="");
  pc.data   <- getURL(pc.url);
  
  pc.data   <- xmlToList(xmlParse(pc.data));
  pc.data   <- as.matrix(pc.data$IdList);
  colnames(pc.data) <- 'CID';
  rownames(pc.data) <- NULL;
  
  return(pc.data);
}

db.pubchem.getIdByFormula <- function(form)
{  
  pc.loc    <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?";
  pc.query  <- "db=pccompound&term=";
  pc.url    <- paste(pc.loc, pc.query, form, sep="");
  pc.data   <- getURL(pc.url);
  
  pc.data   <- xmlToList(xmlParse(pc.data));
  pc.data   <- as.matrix(pc.data$IdList);
  colnames(pc.data) <- 'CID';
  rownames(pc.data) <- NULL;
  
  return(pc.data);
}

db.pubchem.getMolecule <- function(pc.id)
{
  pc.loc    <- "http://pubchem.ncbi.nlm.nih.gov/rest/pug/";
  pc.mol    <- NULL;
  
  for (i in (1:length(pc.id)))
  {
    pc.input  <- paste("compound/cid/",pc.id[i,],"/SDF",sep="");
    pc.url    <- paste(pc.loc,pc.input,sep="");
    pc.mol    <- c(pc.mol, load.molecules(molfiles=pc.url));    
  }
  
  return(pc.mol);
}