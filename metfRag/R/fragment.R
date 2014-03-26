frag.loadFragmenter <- function(frag.path="C:/Users/Adrian Helmchen/Documents/MetFrag/lib/")
{
  .jaddClassPath(frag.path);
  fobj <- new(J("de/ipbhalle/metfrag/fragmenter/Fragmenter"), F, T, F);
  
  return(fobj);
}

frag.smiles <- function(smiles=NULL)
{
  m <- parse.smiles("CN(C)CC(C1=C=C(C=C1)OC)C2(CCCCC2)O");
  m <- m[[1]];

  J("org/openscience/cdk/tools/manipulator/AtomContainerManipulator")$percieveAtomTypesAndConfigureAtoms(m);
  convert.implicit.to.explicit(m);

  return(m);
}

frag.generateFragments <- function(path=NULL, smiles=NULL)
{
  fobj <- frag.loadFragmenter(path);
  smiles <- frag.smiles(smiles);
  
  v <- new(J("java/util/Vector"));
  fobj$setPeakList(v);
  
  frags <- fobj$generateFragmentsInMemory(smiles, F, as.integer(3));
  frags <- as.list(frags);
  
  return(frags);
}