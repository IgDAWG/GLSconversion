#' Check Input Parameterss
#'
#' Check input parameters for invalid entries.
#' @param Convert String Direction for conversion.
#' @param Output String Type of output.
#' @param System String Genetic system (HLA or KIR) of the data being converted
#' @param HZY.Red Logical Reduction of homozygote genotypes to single allele.
#' @param DRB345.Check Logical Check DR haplotypes for consistency and flag unusual haplotypes.
#' @param Cores.Lim Integer How many cores can be used.
#' @note This function is for internal use only.
Check.Params <- function (Convert,Ouput,System,HZY.Red,DRB345.Check,Cores.Lim) {

  if( is.na(match(Convert,c("GL2Tab","Tab2GL"))) ) { Err.Log("P.Convert") ; stop("Conversion Stopped.",call.=FALSE) }
  if( is.na(match(Output,c("R","txt","csv","pypop"))) ) { Err.Log("P.Output") ; stop("Conversion Stopped.",call.=FALSE) }
  if( is.na(match(System,c("HLA","KIR"))) ) { Err.Log("P.System") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(HZY.Red) ) { Err.Log("P.HZY") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.logical(DRB345.Check) ) { Err.Log("P.DRB") ; stop("Conversion Stopped.",call.=FALSE) }
  if( !is.numeric(Cores.Lim) || !is.integer(Cores.Lim) ) { Err.Log("P.Cores") ; stop("Conversion Stopped.",call.=FALSE) }

}

#' GL String Locus Check
#'
#' Check GL string for locus appear in multiple gene fields.
#' @param x GL String to check against
#' @param Loci Loci to check
#' @note This function is for internal use only.
CheckString <- function(x,Loci) {

  test <- sapply(Loci,FUN = function(z) regexpr(z,x) )
  test[test==-1] <- NA
  test.CS <- colSums(test, na.rm=TRUE)
  if( max(test.CS)>1 ) {

    Loci.Err <- paste(Loci[which(test.CS>1)],collapse=",")
    GLS <- paste(x,collapse="^")
    Err.Log("Locus.MultiField",GLS,Loci.Err)
    stop("Conversion Stopped.",call.=FALSE)

  }

}

