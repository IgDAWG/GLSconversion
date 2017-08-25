#' File Name Extraction
#'
#' Function to extract file path.
#' @param x File name.
#' @note This function is for internal use only.
getName <- function(x) {

  tmpDir <- dirname(x)
  tmpName <- basename(x)

  if(basename(x)==x) {
    outName <- paste("Converted_",x,sep="")
  } else {
    outName <- paste(tmpDir,"/Converted_",tmpName,sep="")
  }

  return(outName)

}

#' Build Output Matrix for GL2Tab Conversion
#'
#' Initializes output matrix format for GL2Tab conversion
#' @param System Character Genetic system HLA- or KIR
#' @param Loci The loci for header names
#' @note This function is for internal use only.
Build.Matrix <- function(System,Loci) {

  if( sum(grepl("DRB.HapFlag",Loci))>0 ) { Loci <- Loci[-grep("DRB.HapFlag",Loci)] }
  Loci.Grp <- rep(Loci,each=2)

  if(System=="HLA-") {
    Out <- mat.or.vec(nr=1,nc=length(Loci.Grp)+1) ; colnames(Out) <- c(Loci.Grp,"DRB.HapFlag")
  } else {
    Out <- mat.or.vec(nr=1,nc=length(Loci.Grp)) ; colnames(Out) <- Loci.Grp
  }
  colnames(Out)[seq(1,length(Loci.Grp),by=2)] <- paste(Loci,"_1",sep="")
  colnames(Out)[seq(2,length(Loci.Grp),by=2)] <- paste(Loci,"_2",sep="")

  return(Out)

}

#' Tabular Data Locus Format Tool
#'
#' Correctly orders the expanded GL string
#' @param x Single row of converted GL string
#' @param Tab.Out Single row data frame for mapping converted GL strings
#' @note This function is for internal use only.
Format.Tab <- function(x,Tab.Out) {

  Tab.Out[,match(colnames(x),colnames(Tab.Out))] <- x
  return(Tab.Out)

}

#' Remove Locus Names for Ambiguous Alleles
#'
#' Remove Locus name for each allele in an ambiguous allele string
#' @param x Allele
#' @note This function is for internal use only.
Format.Allele <- function(x) {

  if(grepl("/",x)) {
    tmp <- strsplit(unlist(strsplit(x,"/")),"\\*")
    tmp <- paste(unlist(lapply(tmp,"[",1)[1]),
           paste(unlist(lapply(tmp,"[",2)),collapse="/"),
           sep="*")
  } else { tmp <- x }
  return(tmp)

}
