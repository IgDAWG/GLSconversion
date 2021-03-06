#' DRB345 haplotype zygosity checker
#'
#' Checks DR haplotypes for correct zygosity and flags unanticipated haplotypes
#' @param Locus Locus of interest to test for consistency
#' @param Genotype Row of data set data frame following DRB345 parsing
#' @note This function is for internal use only.
DRB345.Check.Zygosity <- function(Locus,Genotype) {

  Genotype <- Filler(Genotype,Type="Remove") ; Genotype <- Genotype[which(Genotype!="")]

  DR.out <- data.frame(Locus_1=character(), Locus_2=character(), Flag=character(), stringsAsFactors=F)
  Abs <- paste(Locus,"*00:00",sep="")

  DR.Locus <- gsub("HLA-","",Locus) ; DR.Calls <- gsub("HLA-","",Genotype)
  DR.Calls <- sapply(DR.Calls,FUN=GetField.GLS,Res=1) # get 1 Field Resolution for Genotype Calls
  names(DR.Calls) <- NULL ; Flag <- NULL

  #DRB1 - get expected DRB3/4/5 genotypes
  DR345.Exp.Calls <- DRB345.Exp(DR.Calls[grep("DRB1",DR.Calls)])

  #DRB345 Check
  getDRB345 <- grep(DR.Locus,DR.Calls) ; DR.obs <- length(getDRB345) ; DR.exp <- sum(grepl(DR.Locus,DR345.Exp.Calls))

  # Assign Genotypes
  if( DR.obs != DR.exp ) {

    # Inconsistent Genotype Possibilities
    if ( DR.obs==0 && DR.exp>=1 ) {
      # Missing Allele
      DR.out[1, 'Locus_1'] <- Abs ; DR.out[1, 'Locus_2'] <- Abs ; DR.out[1, 'Flag'] <- paste(Locus,"_M",sep="")
    } else if ( DR.obs >=1 && DR.exp==0 ) {
      # Extra Allele
      DR.out[1, 'Locus_1'] <- Genotype[getDRB345[1]] ; DR.out[1, 'Locus_2'] <- Abs ; DR.out[1, 'Flag'] <- paste(Locus,"_P",sep="")
    } else if( DR.obs==1 && DR.exp==2 ) {
      # Presumed Homozygote Missing Allele
      DR.out[1, 'Locus_1'] <- Genotype[getDRB345[1]] ; DR.out[1, 'Locus_2'] <- Genotype[getDRB345[1]] ; DR.out[1, 'Flag'] <- ""
    } else if( DR.obs==2 && DR.exp==1 ) {

      if( Genotype[getDRB345[1]] == Genotype[getDRB345[2]] ) {
        # Extra Allele ... False Homozygote Assumption
        DR.out[1, 'Locus_1'] <- Genotype[getDRB345[1]] ; DR.out[1, 'Locus_2'] <- Abs ; DR.out[1, 'Flag'] <- ""
      } else {
        # Extra Allele Present
        DR.out[1, 'Locus_1'] <- Genotype[getDRB345[1]] ; DR.out[1, 'Locus_2'] <-Genotype[getDRB345[2]] ; DR.out[1, 'Flag'] <- paste(Locus,"_P",sep="")
      }

    }

  } else {

    DR.out[1, 'Flag'] <- ""

    # Consistent Genotype
    if(  DR.obs==0 ) {
      DR.out[1, 'Locus_1'] <-Abs ; DR.out[1, 'Locus_2'] <- Abs
    } else if( DR.obs==1 ) {
      DR.out[1, 'Locus_1'] <- Genotype[getDRB345[1]] ; DR.out[1, 'Locus_2'] <- Abs
    } else if ( DR.obs==2 ) {
      DR.out[1, 'Locus_1'] <- Genotype[getDRB345[1]] ; DR.out[1, 'Locus_2'] <- Genotype[getDRB345[2]]
    }

  }

  # Return Result
  return(DR.out)

}

#' DRB345 Expected
#'
#' Checks DRB1 Genotype and Returns Expected DR345 Loci
#' @param DRB1.Genotype DRB1 Subject Genotypes
#' @note This function is for internal use only.
DRB345.Exp <- function(DRB1.Genotype) {

  #Checks for and fixes certain DRB345 errors that are consistent with known DR haplotypes
  Rules <- list("DRB1*01"="","DRB1*10"="","DRB1*08"="",
                "DRB1*03"="DRB3","DRB1*11"="DRB3","DRB1*12"="DRB3","DRB1*13"="DRB3","DRB1*14"="DRB3",
                "DRB1*04"="DRB4","DRB1*07"="DRB4","DRB1*09"="DRB4",
                "DRB1*15"="DRB5","DRB1*16"="DRB5")

  DRB1.Genotype <- gsub("HLA-","",DRB1.Genotype)
  DRB1.Genotype <- sapply(DRB1.Genotype,FUN=GetField.GLS,Res=1)

  # Allele 1
  DRB1.1 <- DRB1.Genotype[1]
  DR.Gtype <- as.character(Rules[DRB1.1])

  # Allele 2
  if(length(DRB1.Genotype)==1) {
    #Consider Homozygote
    DRB1.2 <-  DRB1.Genotype[1]
  } else {
    DRB1.2 <-  DRB1.Genotype[2]
  }
  DR.Gtype <- c(DR.Gtype,as.character(Rules[DRB1.2]))
  DR.Gtype <- DR.Gtype[which(DR.Gtype!="")]

  if(length(DR.Gtype)>0) { DR.Gtype <- paste("HLA-",DR.Gtype,sep="") }
  return(DR.Gtype)

}

#' HLA trimming function
#'
#' Trim a properly formatted unambiguous HLA allele to desired number of fields.
#' @param x HLA allele.
#' @param Res Resolution desired.
#' @note This function is for internal use only.
GetField.GLS <- function(x,Res) {
  Tmp <- unlist(strsplit(as.character(x),":"))
  if (length(Tmp)<2) {
    return(x)
  } else if (Res==1) {
    return(Tmp[1])
  } else if (Res > 1) {
    Out <- paste(Tmp[1:Res],collapse=":")
    return(Out)
  }
}
