#' DRB345 haplotype zygosity checker
#'
#' Checks DR haplotypes for correct zygosity and flags unanticipated haplotypes
#' @param Locus Locus of interest to test for consistency
#' @param Genotype Row of data set data frame following DRB345 parsing
#' @note This function is for internal BIGDAWG use only.
DRB345.Check.Zygosity <- function(Locus,Genotype) {

  #Checks for and fixes certain DRB345 errors that are consistent with known DR haplotypes
  Rules <- list("DRB1*01"="^","DRB1*10"="^","DRB1*08"="^",
                "DRB1*03"="DRB3","DRB1*11"="DRB3","DRB1*12"="DRB3","DRB1*13"="DRB3","DRB1*14"="DRB3",
                "DRB1*04"="DRB4","DRB1*07"="DRB4","DRB1*09"="DRB4",
                "DRB1*15"="DRB5","DRB1*16"="DRB5")

  DR.out <- data.frame(Locus_1=character(), Locus_2=character(), Flag=character(), stringsAsFactors=F)
  Abs <- paste(Locus,"*^",sep="")

  DR.Calls <- gsub("HLA-","",Genotype) ; DR.Locus <- gsub("HLA-","",Locus)
  DR.Calls <- sapply(DR.Calls,FUN=GetField,Res=1) # get 1 Field Resolution for Genotype Calls
  names(DR.Calls) <- NULL ; Flag <- NULL

  #DRB1 - get expected DRB3/4/5 genotypes
  getDRB1 <- grep("DRB1",DR.Calls)

  DRB1.1 <- DR.Calls[getDRB1[1]]
  DR.Gtype <- as.character(Rules[DRB1.1])

  if(length(getDRB1)==1) {
    DRB1.2 <-  DR.Calls[getDRB1[1]]
  } else {
    DRB1.2 <-  DR.Calls[getDRB1[2]]
  }
  DR.Gtype <- c(DR.Gtype,as.character(Rules[DRB1.2]))

  #DRB345 Check
  getDRB345 <- grep(DR.Locus,DR.Calls)
  if( length(getDRB345)>0 ) {

    DR.obs <- length(getDRB345)
    DR.exp <- sum(grepl(DR.Locus,DR.Gtype))
    A1 <- Genotype[getDRB345[1]] ; A2 <- ifelse(length(getDRB345)>1, Genotype[getDRB345[2]], Abs)

    if( DR.obs != DR.exp ) {
      if( DR.obs==1 && DR.exp==2 ) {
        DR.out[1, 'Locus_1'] <- A1 ; DR.out[1, 'Locus_2'] <- A1 ; DR.out[1, 'Flag'] <- F
      } else if( DR.obs==2 && DR.exp==1 && A1==A2 ) {
        DR.out[1, 'Locus_1'] <- A1 ; DR.out[1, 'Locus_2'] <- Abs ; DR.out[1, 'Flag'] <- F
      } else if( DR.obs==2 && DR.exp==1 && A1!=A2 ) {
        DR.out[1, 'Locus_1'] <- A1 ; DR.out[1, 'Locus_2'] <- A2 ; DR.out[1, 'Flag'] <- T
      }
    } else {
      DR.out[1, 'Locus_1'] <- A1 ; DR.out[1, 'Locus_2'] <- A2 ; DR.out[1, 'Flag'] <- F
    }

  } else { DR.out[1, 'Locus_1'] <- Abs ; DR.out[1, 'Locus_2'] <- Abs ; DR.out[1, 'Flag'] <- F }

  # Return Result
  return(DR.out)

}

#' HLA trimming function
#'
#' Trim a properly formatted HLA allele to desired number of fields.
#' @param x HLA allele.
#' @param Res Resolution desired.
#' @note This function is for internal BIGDAWG use only.
GetField <- function(x,Res) {
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
