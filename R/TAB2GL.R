#' Genotype List String to Tabular Data Conversion
#'
#' Expands GL strings to columns of adjacent locus pairs.
#' @param df Data frame containing GL strings
#' @param System Character Genetic system HLA or KIR
#' @param HZY.Red Logical Should homozygote genotypes be a single allele for non-DRB345.
#' @param Cores Integer How many cores can be used.
Tab2GL.conv <- function(df,System,HZY.Red,Cores) {

  # Check for ambiguous data at allele ("/")
  if( sum(grepl("\\|",df[,3]))>0 ) { stop("This appears to be ambiguous data. Conversion stopped.",call.=F) }

  # Check for column formatting consistency
  if( ncol(df) < 4 ) { stop("Your data is not properly formatted for the Tab2GL parameter. Conversion stopped.",call.=F) }
  if( !is.even(ncol(df)) )  { stop("Your data is not properly formatted for the Tab2GL parameter. Conversion stopped.",call.=F) }

  # Run Conversion
  df.list <- lapply(seq(1,nrow(df)),FUN= function(i) df[i,3:ncol(df)])
  GL <- mclapply(df.list,FUN=Tab2GL,System=System,HZY.Red=HZY.Red,mc.cores=Cores)
  GL <- do.call(rbind,GL)
  colnames(GL) <- "GL.String"
  GL <- cbind(df[,1:2],GL)

  return(GL)

}

#' Genotype List String Condenser
#'
#' Condenses column of loci into a GL string
#' @param x Row of loci to condense
#' @param System Character Genetic system HLA or KIR
#' @param HZY.Red Logical Should homozygote genotypes be a single allele for non-DRB345.
Tab2GL <- function(x,System,HZY.Red) {

  x <- x[which(x!="")]
  Loci <- unique(sapply(colnames(x),FUN=gsub,pattern="\\.1|\\.2|\\_1|\\_2",replacement=""))

  GLS <- NULL
  # Condense Alleles (+)
  for(i in Loci) {

    Alleles <- x[grep(i,colnames(x))]
    if( sum(grepl(System,Alleles))==0 ) { Alleles <- paste(System,i,"*",Alleles,sep="") }
    A1 <- Alleles[1] ; A2 <- Alleles[2]
    if( is.na(A1) || is.na(A2) ) {
      GLS <- c(GLS,Alleles)
    } else if( HZY.Red && A1==A2 ) { GLS <- c(GLS,Alleles[1])
    } else { GLS <- c(GLS,paste(Alleles,collapse="+"))  }

  }

  # Condense Chromosomes (^)
  GLS <- paste(GLS,collapse="^")

  return(GLS)

}

