#' Genotype List String to Tabular Data Conversion
#'
#' Expands GL strings to columns of adjacent locus pairs.
#' @param df Data frame containing GL strings
#' @param System Character Genetic system HLA or KIR
#' @param Cores Integer How many cores can be used.
#' @note This function is for internal use only.
GL2Tab.wrapper <- function(df,System,Cores) {

  # Data column
  LastCol <- ncol(df)
  MiscCol <- seq(1,ncol(df)-1)

  # Run Conversion
  df.list <- strsplit(df[,LastCol],"\\^")
  Tab <- parallel::mclapply(df.list,FUN=GL2Tab,System=System,mc.cores=Cores)
    Loci <- sort(unique(gsub("_1|_2","",unlist(lapply(Tab,colnames)))))
    Tab.Out <- Build.Matrix(System,Loci)

  Tab <- parallel::mclapply(Tab,FUN=Format.Tab,Tab.Out=Tab.Out,mc.cores=Cores)
    Tab <- do.call(rbind,Tab)
    Tab[Tab==0] <- ""
    Tab[grepl("\\^",Tab)] <- ""

  Tab <- cbind(df[,MiscCol],Tab)
  return(Tab)

}

#' Genotype List String Expander
#'
#' Expands GL string into a table of adjacent loci
#' @param x Character GL string to expand
#' @param System Character Genetic system HLA or KIR
#' @note This function is for internal use only.
GL2Tab <- function(x,System) {

  # Break GL String and Remove Any Absent Call Type Strings (00:00)
  Calls <- unlist(sapply(x,FUN=function(x) strsplit(x,"\\+"))) ; names(Calls) <- NULL

  # Check GL String For Locus*Allele/Locus*Allele Ambiguity Formatting
  invisible(sapply(Calls,CheckString.Allele))

  # Collapse Ambiguous allele names to remove locus prefix
  if( sum(grepl("/",Calls)>0 ) ) {
    Calls <- unlist(lapply(Calls,Format.Allele,Type="off"))
  }

  # Get loci and initialize table
  Loci <- unique(unlist(lapply(strsplit(Calls,"\\*"),"[",1)))
  if(System=="HLA-") {
    if( sum(grepl("DRB1",Loci))>0 ) { Loci <- c(Loci,DRB345.Exp(Calls[grep("DRB1",Calls)])) }
    if(sum(grepl("\\^",Loci))>0) { Loci <- Loci[-grep("\\^",Loci)] }
    Loci <- unique(Loci)
  }
  Tab <- Build.Matrix(System,Loci)

  # Check GL String For Locus^Gene Field Consistency
  invisible(CheckString.Locus(x,Loci))

  # Populate table
  if(System=="HLA-") { DRB345.Flag <- NULL }
  for(i in Loci) {

    if(System=="HLA-") {
      # HLA System
      if(i=="HLA-DRB3" || i=="HLA-DRB4" || i=="HLA-DRB5") {
        if( sum(grepl("DRB1",x))>0 ) {
          # Assumptions for DRB345
          DRB.GTYPE <- DRB345.Check.Zygosity(i,Calls[grep("DRB",Calls)])
          Tab[1,grep(i,colnames(Tab))] <- as.character(DRB.GTYPE[1,c('Locus_1','Locus_2')])
          # for inconsistent DR haplotypes
          if( as.logical(DRB.GTYPE[,'Flag']) ) { DRB345.Flag <- c(DRB345.Flag,i) }
        } else {
          # No DRB1 but DBR345 (ZYgosity Check Not Determined)
          if( sum(grepl("DRB",x))>0 ) { DRB345.Flag <- "ND" }
          Tab[1,grep(i,colnames(Tab))] <- Calls[grep(i,Calls)]
        }
      } else {
        # no DRB345
        Tab[1,grep(i,colnames(Tab))] <- Calls[grep(i,Calls)]
      } # fi DRB345
    } else {
      # non-HLA System
      Tab[1,grep(i,colnames(Tab))] <- Calls[grep(i,Calls)]
    }

  if(System=="HLA-") { Tab[1,'DRB.HapFlag'] <- ifelse(!is.null(DRB345.Flag), paste(unlist(DRB345.Flag),collapse=",") , "") }

  }# fi Loci Loop

  return(Tab)

}

