#' Genotype List String to Tabular Data Conversion
#'
#' Expands GL strings to columns of adjacent locus pairs.
#' @param df Data frame containing GL strings
#' @param System Character Genetic system HLA or KIR
#' @param HZY.Red Logical Should homozygote genotypes be a single allele for non-DRB345.
#' @param Cores Integer How many cores can be used.
#' @note This function is for internal use only.
Tab2GL.wrapper <- function(df,System,HZY.Red,Cores) {

  # Define data locus columns assuming Locus columns come in pairs
  colnames(df) <- sapply(colnames(df),FUN=gsub,pattern="\\.1|\\.2|\\_1|\\_2",replacement="")
  DataCol <- as.numeric(sapply(names(which(table(colnames(df))==2)), FUN=function(x) grep(x,colnames(df))))
  MiscCol <- setdiff(1:ncol(df),DataCol)

  # Check for identical rows of miscellanous information from non-data columns. Ambiguous Data Flag.
  Misc.tmp <- apply(df[,MiscCol],MARGIN=1,FUN=paste,collapse=":")
  if( length(which(table(Misc.tmp)>1))>0 ) {
    Err.Log("Table.Amb") ; stop("Conversion stopped.",call.=F)
  }; rm(Misc.tmp)

  # Pre-format data to SystemLoci*Allele if necessary
  if( sum(grepl(System,colnames(df)[DataCol]))==0 ) { colnames(df)[DataCol] <- paste(System,colnames(df)[DataCol],sep="") }
  for(i in DataCol) {
    if(sum(grepl(colnames(df)[i],df[,i]))==0) {
      df[,i] <- sapply(df[,i], FUN = Append.System, df.name=colnames(df)[i] )
    }
  }

  # Run Conversion
  df.list <- lapply(seq(1,nrow(df)), FUN=function(k) df[k,DataCol])
  GL <- parallel::mclapply(df.list,FUN=Tab2GL,System=System,HZY.Red=HZY.Red,mc.cores=Cores)
  GL <- do.call(rbind,GL)

  if(ncol(GL)==1) { colnames(GL) <- "GL.String" } else if(ncol(GL)==2) { colnames(GL) <- c("GL.String","DRB.HapFlag") }
  GL <- cbind(df[,MiscCol],GL)

  return(GL)

}

#' Genotype List String Condenser
#'
#' Condenses column of loci into a GL string
#' @param x Row of loci to condense
#' @param System Character Genetic system HLA or KIR
#' @param HZY.Red Logical Should homozygote genotypes be a single allele for non-DRB345.
#' @note This function is for internal use only.
Tab2GL <- function(x,System,HZY.Red) {

  #x <- rmABstrings(x)
  x <- x[which(x!="")]
  colnames(x) <- sapply(colnames(x),FUN=gsub,pattern="\\.1|\\.2|\\_1|\\_2",replacement="")
  Loci <- unique(colnames(x))
  if(System=="HLA-") {
    if( sum(grepl("DRB1",Loci))>0 ) { Loci <- c(Loci,DRB345.Exp(x[grep("DRB1",colnames(x))])) }
    if( sum(grepl("\\^",Loci))>0 ) { Loci <- Loci[-grep("\\^",Loci)] }
    Loci <- unique(Loci)
  }

  # Append Locus to ambiguous Allele/Allele calls
  x[] <- sapply(x,Format.Allele,Type="on")

  GLS <- NULL ; if(System=="HLA-") { DRB345.Flag <- NULL }
  # Condense Alleles (+)
  for(i in Loci) {

    Alleles <- as.character(x[,grep(i,colnames(x))])
    A1 <- Alleles[1] ; A2 <- Alleles[2]

    if(System=="HLA-") {
      if(i=="HLA-DRB3" || i=="HLA-DRB4" || i=="HLA-DRB5") {
          if( sum(grepl("DRB1",x))>0 ) {
            DRB.GTYPE <- DRB345.Check.Zygosity(i, x[grep("DRB",x)] )
            DRB.GTYPE[grepl("\\^",DRB.GTYPE)] <- NA
            A1 <- DRB.GTYPE[,'Locus_1'] ; A2 <- DRB.GTYPE[,'Locus_2'] ; Alleles <- c(A1,A2)
            if(DRB.GTYPE$Flag) {
              # DRB345 is not consistent
              if( as.logical(DRB.GTYPE[,'Flag']) ) { DRB345.Flag <- c(DRB345.Flag,i) }
            }
            # DRB345 but no DRB1 (ZYgosity Check Not Determined)
          } else if( sum(grepl(grep("DRB",x)))>0 ) { DRB345.Flag <- "ND" }
      } # fi DRB345
    } # fi HLA

    if( is.na(A1) || is.na(A2) ) {
      GLS <- c(GLS,Alleles)
    } else if( HZY.Red && A1==A2 ) { GLS <- c(GLS,Alleles[1])
    } else { GLS <- c(GLS,paste(Alleles,collapse="+"))  }

  } # fi Loci Loop

  # Condense Chromosomes (^)
  GLS <- paste(GLS,collapse="^")

  if(System=="HLA-") {
    DRB.HapFlag <- ifelse(!is.null(DRB345.Flag), paste(unlist(DRB345.Flag),collapse=",") , "")
    Out <- c(GLS,DRB.HapFlag)
  } else {
    Out <- GLS
  }

  return(Out)

}
