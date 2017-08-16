#' Genotype List String to Tabular Data Conversion
#'
#' Expands GL strings to columns of adjacent locus pairs.
#' @param df Data frame containing GL strings
#' @param System Character Genetic system HLA or KIR
#' @param DRB345.Check Logical Flag unusual DR haplotypes.
#' @param Cores Integer How many cores can be used.
GL2Tab.conv <- function(df,System,DRB345.Check,Cores) {

  # Check for ambiguous data at Locus "/" or genotype "|"
  if( sum(grepl("\\|",df[,3]))>0 ) { stop("This appears to be ambiguous data. Conversion stopped.",call.=F) }
  if( sum(grepl("\\/",df[,3]))>0 ) { stop("This appears to be ambiguous data. Conversion stopped.",call.=F) }

  # Run Conversion
  Tab <- mclapply(df[,3],FUN=GL2Tab,System=System,DRB345.Check=DRB345.Check,mc.cores=Cores)
    Loci <- sort(unique(gsub("_1|_2","",unlist(lapply(Tab,colnames)))))
    Loci.Grp <- rep(Loci,each=2)
    Out <- mat.or.vec(nr=1,nc=length(Loci.Grp)) ; colnames(Out) <- Loci.Grp
    colnames(Out)[seq(1,length(Loci.Grp),by=2)] <- paste(Loci,"_1",sep="")
    colnames(Out)[seq(2,length(Loci.Grp),by=2)] <- paste(Loci,"_2",sep="")

  Tab <- mclapply(Tab,FUN=format.Tab,Out=Out,mc.cores=Cores)
    Tab <- do.call(rbind,Tab)
    Tab[Tab==0] <- ""
    Tab[grepl("\\^",Tab)] <- ""

  Tab <- cbind(df[,c(1,2)],Tab)
  return(Tab)

}

#' Genotype List String Expander
#'
#' Expands GL string into a table of adjacent loci
#' @param x Character GL string to expand
#' @param System Character Genetic system HLA or KIR
#' @param DRB345.Check Logical Flag unusual DR haplotypes.
GL2Tab <- function(x,System,DRB345.Check) {

  # Break GL String
  tmp <- unlist(strsplit(x,"\\^")) # Locus
  tmp <- sapply(tmp,FUN=function(x)strsplit(x,"\\+")) # Chromosome
  Calls <- unlist(tmp) ; names(Calls) <- NULL

  # Get Loci and Initialize Table
  Loci <- unique(unlist(lapply(strsplit(Calls,"\\*"),"[",1)))
  Loci.Grp <- rep(Loci,each=2)
  Tab <- mat.or.vec(nr=1,nc=length(Loci.Grp)) ; colnames(Tab) <- Loci.Grp
  colnames(Tab)[seq(1,length(Loci.Grp),by=2)] <- paste(Loci,"_1",sep="")
  colnames(Tab)[seq(2,length(Loci.Grp),by=2)] <- paste(Loci,"_2",sep="")

  # Populate Table
  for(i in Loci) {

    getCalls <- grep(i,Calls)

    if(System=="HLA-") {
      # Assumptions for DRB345
      if(i=="HLA-DRB3" || i=="HLA-DRB4" || i=="HLA-DRB5") {

        DRB.GTYPE <- DRB345.zygosity(i,Calls[grep("DRB",Calls)])
        if(DRB.GTYPE[,'Flag']) {
          # DRB345 is not consistent
          if(DRB345.check) { Tab[1,grep(i,colnames(Tab))] <- paste(as.character(DRB.GTYPE[1,c('Locus_1','Locus_2')]),"!",sep="")
          } else { Tab[1,grep(i,colnames(Tab))] <- as.character(DRB.GTYPE[1,c('Locus_1','Locus_2')])}
        } else {
          # DBR345 is consistent
          Tab[1,grep(i,colnames(Tab))] <- as.character(DRB.GTYPE[1,c('Locus_1','Locus_2')])
        }
      } else {
        # non-DRB345 situations
        Tab[1,grep(i,colnames(Tab))] <- Calls[grep(i,Calls)]
      }

    } else {
      # non-HLA
      Tab[1,grep(i,colnames(Tab))] <- Calls[grep(i,Calls)]
    }

  }# End i

  return(Tab)

}

#' Tabular Data Locus Format Tool
#'
#' Correctly orders the expanded GL string
#' @param x Single row of converted GL string
#' @param Out Single row data frame for mapping converted GL strings
#' @note This function is for internal BIGDAWG use only.
format.Tab <- function(x,Out) {

  Out[,match(colnames(x),colnames(Out))] <- x
  return(Out)

}
