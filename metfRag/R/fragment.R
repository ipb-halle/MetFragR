.packageName <- "metfRag"

require(rJava, quietly=TRUE)
require(rcdk, quietly=TRUE)

.onLoad<-function(libname, pkgname) {
	jar.metfrag <- paste(libname, pkgname, "java", "MetFragR-2.4.2-jar-with-dependencies.jar", sep=.Platform$file.sep)	
	.jinit(classpath=c(jar.metfrag))
}

frag.generateFragments <- function(molecule, treeDepth = 2)
{
  if (missing(molecule) == TRUE)
  {
    stop("ERROR: Molecule is missing.");
  }

  if(class(molecule) == "character") container=parse.smiles(molecule)
  else container=molecule

  obj = .jnew("de/ipbhalle/metfrag/r/MetfRag")
  frags<- .jcall(obj, '[Lorg/openscience/cdk/interfaces/IAtomContainer;',
               'generateAllFragments', 
               container, as.integer(treeDepth))

  frags <- as.list(frags);
  return(frags);
}

frag.generateMatchingFragments <- function(molecule, mzs, neutralMonoisotopicMass, mzabs = 0.01, mzppm = 10.0, posCharge = TRUE, ionMode = 1, treeDepth = 2)
{
  if (missing(molecule) == TRUE)
  {
    stop("ERROR: Molecule is missing.");
  }
  if (missing(mzs) == TRUE)
  {
    stop("ERROR: MZ values are missing.");
  }

  if(class(molecule) == "character") container=parse.smiles(molecule) 
  else container=molecule

  obj = .jnew("de/ipbhalle/metfrag/r/MetfRag")
  frags<- .jcall(obj, '[Lorg/openscience/cdk/interfaces/IAtomContainer;',
               'generateMatchingFragments',         
               container, .jarray(mzs,"[D"), neutralMonoisotopicMass, as.double(mzabs), as.double(mzppm), as.logical(posCharge), as.integer(ionMode), as.integer(treeDepth));

  frags <- as.list(frags);
  return(frags);
}
