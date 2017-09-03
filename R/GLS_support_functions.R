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

  outName <- gsub(".txt","",outName)
  return(outName)

}

#' Replace 00:00 allele strings
#'
#' Replaces 00:00 absent allele strings with a blank.
#' @param x Genotype Column.
#' @param Locus Locus column to adjust.
#' @param Type String specifying whether to pad ('Fill') or leave blank ('Remove') absent calls
#' @note This function is for internal use only.
Filler <- function(x,Locus=NULL,Type) {

  if (Type=="Fill") {
    which(x=="")
    Locus <- gsub("_1|_2","",Locus)
    x[which(x=="")] <- paste(Locus,"*00:00",sep="")

  }

  if (Type=="Remove") {

    x[grep("\\*00",x)] <- ""

  }

  return(x)

}

#' Build Output Matrix for GL2Tab Conversion
#'
#' Initializes output matrix format for GL2Tab conversion
#' @param System Character Genetic system HLA- or KIR
#' @param Loci The loci for header names
#' @note This function is for internal use only.
Build.Matrix <- function(System,Loci) {

  Loci.Grp <- rep(Loci,each=2)

  if(System=="HLA-") {
    Out <- mat.or.vec(nr=1,nc=length(Loci.Grp)+1) ; colnames(Out) <- c(Loci.Grp,"DR.HapFlag")
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
#' @param Order Single row data frame for mapping converted GL strings
#' @note This function is for internal use only.
Format.Tab <- function(x,Order) {

  Order[,match(colnames(x),colnames(Order))] <- x
  return(Order)

}

#' Remove or Append Locus Names for Ambiguous Alleles
#'
#' Remove or Append Locus name from/to allele in an ambiguous allele string
#' @param x Allele String
#' @param Type String specifying whether to strip ('off') or append ('on') locus prefix
#' @note This function is for internal use only.
Format.Allele <- function(x,Type) {

  if(Type=="off") {
    if(grepl("/",x)) {
      tmp <- strsplit(unlist(strsplit(x,"/")),"\\*")
      Fix <- paste(unlist(lapply(tmp,"[",1)[1]),
             paste(unlist(lapply(tmp,"[",2)),collapse="/"),
             sep="*")
    } else {  Fix <- x }
  }

  if(Type=="on"){
   if(grepl("/",x)) {
      Locus <- unlist(strsplit(x,"\\*"))[1]
      Fix <- paste(
              paste(
               Locus,unlist(strsplit(unlist(strsplit(x,"\\*"))[2],"/"))
               ,sep="*")
             ,collapse="/")
   } else { Fix <- x }
  }

  return(Fix)
}

#' Append Genetic System Locus Designation to Allele String
#'
#' Adds genetic system (HLA/KIR) to each allele name
#' @param x Vector Column genotypes to append
#' @param df.name String SystemLocus name for each allele.
#' @note This function is for internal use only.
Append.System <- function(x,df.name) {

  getAllele <- which(x!="")
  x[getAllele] <- paste(df.name,x[getAllele],sep="*")
  return(x)

}

#' Removes System and Locus from Alleles
#'
#' Removes the System and Locus designations for alleles calls in GL2Tab
#' @param x Allele
#' @note This function is for internal use only.
Stripper <- function(x) {

  if(x!="") { Fix <- unlist(strsplit(x,"\\*"))[2] } else { Fix <- "" }
  return(Fix)

}
