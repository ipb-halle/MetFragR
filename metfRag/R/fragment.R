frag.loadFragmenter <- function(frag.path)
{  
  .jaddClassPath(frag.path);
  fobj <- new(J("de/ipbhalle/metfrag/fragmenter/Fragmenter"), F, T, F);
  
  return(fobj);
}

frag.smiles <- function(smiles)
{
  m <- parse.smiles(smiles);
  m <- m[[1]];

  J("org/openscience/cdk/tools/manipulator/AtomContainerManipulator")$percieveAtomTypesAndConfigureAtoms(m);
  convert.implicit.to.explicit(m);

  return(m);
}

frag.generateFragments <- function(path, smiles)
{
  if (missing(path) == TRUE)
  {
    #stop("A path to MetFrag binaries have to be given.");
    path <- "D:/Documents/MetFrag/lib/";
  }
  
  if (missing(smiles) == TRUE)
  {
    #stop("A SMILES string for a compound have to be given.");
    smiles <- "CN(C)CC(C1=C=C(C=C1)OC)C2(CCCCC2)O";
  }  
  
  fobj <- frag.loadFragmenter(path);
  smiles <- frag.smiles(smiles);
  
  v <- new(J("java/util/Vector"));
  fobj$setPeakList(v);
  
  frags <- fobj$generateFragmentsInMemory(smiles, F, as.integer(3));
  frags <- as.list(frags);
  
  return(frags);
}