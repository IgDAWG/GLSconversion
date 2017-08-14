#' DRB345 haplotype zygosity checker
#'
#' Checks DR haplotypes for correct zygosity and flags unanticipated haplotypes
#' @param x Row of data set data frame following DRB345 parsing
#' @note This function is for internal BIGDAWG use only.
DRB345.zygosity <- function(DR.Locus,DR.Calls) {

  #Checks for and fixes certain DRB345 errors that are consistent with known DR haplotypes
  Rules <- list("DRB1*01"="^","DRB1*10"="^","DRB1*08"="^",
                "DRB1*03"="DRB3","DRB1*11"="DRB3","DRB1*12"="DRB3","DRB1*13"="DRB3","DRB1*14"="DRB3",
                "DRB1*04"="DRB4","DRB1*07"="DRB4","DRB1*09"="DRB4",
                "DRB1*15"="DRB5","DRB1*16"="DRB5")

  DR.Calls <- sapply(DR.Calls,FUN=GetField,Res=1) # get 1 Field Resolution
  names(DR.Calls) <- NULL ; Flag <- NULL

  #DRB1 - get expected DRB3/4/5 genotypes
  getDRB1 <- grep("DRB1",DR.Calls)
  DRB1.1 <- DR.Calls[getDRB1[1]]
  DR.Gtype <- as.character(Rules[DRB1.1])

  if(length(getDRB1)==1) {
    DRB1.2 <-  DR.Calls[getDRB1[1]]
  } else {
    DRB1.2 <-  DR.Calls[getDRB1[1]]
  }
  DR.Gtype <- c(DR.Gtype,as.character(Rules[DRB1.2]))

  #DRB3 Check
  DRB3.col <- grep("DRB3",colnames(x.1F))
  DRB3.abs <- grep("\\^",x.1F[DRB3.col])
  if( sum(is.na(x.1F[DRB3.col]))==0 ) {

    DRB3.obs <- length(which(sapply(x.1F[DRB3.col],nchar)==7))
    DRB3.exp <- as.numeric(sum(grepl("DRB3",DR.Gtype)))
    A1 <- x.1F[DRB3.col[1]] ; A2 <- x.1F[DRB3.col[2]]

    if( DRB3.obs != DRB3.exp ) {
      if( DRB3.obs==2 && DRB3.exp==1 && A1==A2 ) { x.out[DRB3.col[2]] <- "^"  ; DR3.flag <- F
      } else if( DRB3.obs==2 && DRB3.exp==1 && A1!=A2 ) { DR3.flag <- T
      } else if( DRB3.obs==1 && DRB3.exp==2 ) { x.out[DRB3.col[2]] <- x[DRB3.col[1]]  ; DR3.flag <- F
      } else { DR3.flag <- T }
    } else { DR3.flag <- F }

  } else { DR3.flag <- F }

  #DRB4 Check
  DRB4.col <- grep("DRB4",colnames(x.1F))
  DRB4.abs <- grep("\\^",x.1F[DRB4.col])
  if( sum(is.na(x.1F[DRB4.col]))==0 ) {

    DRB4.obs <- length(which(sapply(x.1F[DRB4.col],nchar)==7))
    DRB4.exp <- as.numeric(sum(grepl("DRB4",DR.Gtype)))
    A1 <- x.1F[DRB4.col[1]] ; A2 <- x.1F[DRB4.col[2]]

    if( DRB4.obs != DRB4.exp ) {
      if( DRB4.obs==2 && DRB4.exp==1 && A1==A2 ) { x.out[DRB4.col[2]] <- "^"  ; DR4.flag <- F
      } else if( DRB4.obs==2 && DRB4.exp==1 && A1!=A2 ) { DR4.flag <- T
      } else if( DRB4.obs==1 && DRB4.exp==2 ) { x.out[DRB4.col[2]] <- x[DRB4.col[1]]  ; DR4.flag <- F
      } else { DR4.flag <- T }
    } else { DR4.flag <- F }

  } else { DR4.flag <- F }

  #DRB5 Check
  DRB5.col <- grep("DRB5",colnames(x.1F))
  DRB5.abs <- grep("\\^",x.1F[DRB5.col])
  if( sum(is.na(x.1F[DRB5.col]))==0 ) {

    DRB5.obs <- length(which(sapply(x.1F[DRB5.col],nchar)==7))
    DRB5.exp <- as.numeric(sum(grepl("DRB5",DR.Gtype)))
    A1 <- x.1F[DRB5.col[1]] ; A2 <- x.1F[DRB5.col[2]]

    if( DRB5.obs != DRB5.exp ) {
      if( DRB5.obs==2 && DRB5.exp==1 && A1==A2 ) { x.out[DRB5.col[2]] <- "^"  ; DR5.flag <- F
      } else if( DRB5.obs==2 && DRB5.exp==1 && A1!=A2 ) { DR5.flag <- T
      } else if( DRB5.obs==1 && DRB5.exp==2 ) { x.out[DRB5.col[2]] <- x[DRB5.col[1]]  ; DR5.flag <- F
      } else { DR5.flag <- T }
    } else { DR5.flag <- F }
  } else { DR5.flag <- F }

  # Set Flag
  if(DR3.flag) { Flag <- c(Flag,"DRB3")  }
  if(DR4.flag) { Flag <- c(Flag,"DRB4")  }
  if(DR5.flag) { Flag <- c(Flag,"DRB5")  }
  if(!is.null(Flag)){ names(Flag) <- "DR_Hap_Error" }

  # Return Result
  Out <- list()
  colnames(x.out) <- colnames(x)
  rownames(x.out) <- NULL
  Out[['Tab']] <- x.out
  Out[['Flag']] <- ifelse(is.null(Flag),"",paste(Flag,collapse=","))
  return(Out)

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
