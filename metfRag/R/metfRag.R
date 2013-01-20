.packageName <- "metfRag"

require(rJava, quietly=TRUE)

.onLoad<-function(libname, pkgname) {
	jar.metfrag <- paste(libname, pkgname, "java", "MetFragSlim.jar", sep=.Platform$file.sep)	
	.jinit(classpath=c(jar.metfrag))
}

score.molecules.from.sdf<-function(sdf, mzs, ints, exact.mass, number.threads=1, mz.abs=0.01, mz.ppm=10, search.ppm=10, pos.charge=TRUE, mode=1, tree.depth=2) {
	
	if(missing(sdf)) stop("Error: SDF is missing!")
	if(missing(mzs)) stop("Error: Vector of mass to charge ratios is missing!")
	if(missing(ints)) stop("Error: Vector of intensities is missing!")
	if(missing(exact.mass)) stop("Error: Neutral exact mass of parent molecule is missing!")

	if(file.access(sdf, 4) != 0) stop("Error: Cannot access SDF! Either it does not exist or is not readable.");
	if(!is.vector(mzs) || !is.numeric(mzs) || any(mzs < 0)) stop("Error: The argument mzs must be a vector of numerals!")
	if(!is.vector(ints) || !is.numeric(ints) || any(ints < 0)) stop("Error: The argument ints must be a vector of positive numerals!")
	if(exact.mass < 0) stop("Error: The argument exact.mass must be a positive numeral!")	
	if(number.threads < 0 || number.threads > 8) stop("Error: The argument number.threads must get values from 0 till 8!")
	if(mz.abs < 0) stop("Error: The argument mz.abs must be positive!")	
	if(mz.ppm < 0) stop("Error: The argument mz.ppm must be positive!")
	if(search.ppm < 0) stop("Error: The argument search.ppm must be positive!")
	if(!is.logical(pos.charge)) stop("Error: The argument pos.charge must be a logical value!")
	if(!is.element(mode, c(-1,0,1))) stop("Error: The argument mode must be a value of -1, 0 or 1")
	if(tree.depth < 1 || tree.depth > 5) stop("Error: The argument tree.depth must get values from 1 till 5!")
	
	obj=.jnew("de/ipbhalle/metfrag/r/MetfRag")
	mols<-.jcall(obj, '[Lorg/openscience/cdk/interfaces/IAtomContainer;', 'scoreMoleculesAgainstSpectrum', sdf, .jarray(mzs,"[D"), .jarray(ints,"[D"), exact.mass, as.integer(number.threads), mz.abs, mz.ppm, search.ppm, pos.charge, as.integer(mode), as.integer(tree.depth))
	if(length(mols) == 0) cat("No results generated!","\n")
	return(mols)
}
