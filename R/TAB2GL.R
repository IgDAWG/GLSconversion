#' Genotype List String to Tabular Data Conversion
#'
#' Expands GL strings to columns of adjacent locus pairs.
#' @param df Data frame containing GL strings
#' @param HZY.Red Logical Should homozygote genotypes be a single allele for non-DRB345.
#' @param Cores Integer How many cores can be used.
Tab2GL.conv <- function(df,HZY.Red,Cores) {

  # Check for ambiguous data at Locus "/"
  if( sum(grepl("\\|",df[,3]))>0 ) { stop("This appears to be ambiguous data. Conversion stopped.",call.=F) }

  # ...
  df.list <- lapply(seq(1,nrow(df)),FUN= function(i) df[i,3:ncol(df)])
  GL <- mclapply(df.list,FUN=Tab2GL,HZY.Red=HZY.Red,mc.cores=Cores)
  GL <- do.call(rbind,GL)
  colnames(GL) <- "GL.String"
  GL <- cbind(df[,1:2],GL)

  return(GL)

}

#' Genotype List String Condenser
#'
#' Condenses column of loci into a GL string
#' @param x Row of loci to condense
#' @param HZY.Red Logical Should homozygote genotypes be a single allele for non-DRB345.
Tab2GL <- function(x,HZY.Red) {

  x <- x[which(x!="")]
  Loci <- unique(sapply(colnames(x),FUN=gsub,pattern="\\.1|\\.2|\\_1|\\_2",replacement=""))

  GLS <- NULL
  # Condense Alleles (+)
  for(i in Loci) {

    Alleles <- x[grep(i,colnames(x))]
    Alleles <- paste("HLA-",i,"*",Alleles,sep="")
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

