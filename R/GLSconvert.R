#' Genotype List String Conversion
#'
#' Main Workhorse wrapper for cross converting columnar table to GL string representaion.
#' @param Data String File name or data frame.
#' @param Convert String Direction for conversion.
#' @param Output String Type of output.
#' @param System String Genetic system (HLA or KIR) of the data being converted
#' @param HZY.Red Logical Reduction of homozygote genotypes to single allele.
#' @param DRB345.Flag Logical Flag unusual DR haplotypes.
#' @param Cores.Lim Integer How many cores can be used.
GLS.Convert <- function(Data,Convert,Output="txt",System="HLA",HZY.Red=FALSE,DRB345.Flag=TRUE,Cores.Lim=1L) {

  # MultiCore Limitations
  if (Cores.Lim!=1L) {
    Cores.Max <- as.integer( floor( parallel::detectCores() * 0.9) )
    if(Sys.info()['sysname']=="Windows" && as.numeric(Cores.Lim)>1) {
      stop("Windows cores max exceeded. Conversion stopped.",call. = F)
    } else if( Cores.Lim > Cores.Max ) { Cores <- Cores.Max
    } else { Cores <- Cores.Lim }
  } else { Cores <- Cores.Lim }

  # Nomenclature system
  if( System == "HLA" ) { System <- "HLA-" }

  # Read in Data and Set Output File Name
  if( is.character(Data) ) {
    if( file.exists(Data) ) {
      df <- read.table(file=Data,header=T,sep="\t",stringsAsFactors=FALSE)
       FP <- getPath(Data)$path
       if( is.na(FP) ) { fileName <- "Converted.txt" } else { fileName <- paste(FP,"Converted.txt",sep="/") }
    } else { stop("Conversion utility cannot local file, please check name. Conversion Stopped.",call.=FALSE) }
  } else { df <- Data ; fileName <- "Converted.txt" }

  # Run Conversion
  switch(Convert,
         GL2Tab = { data.out <- GL2Tab.conv(df,System,DRB345.Flag,Cores) } ,
         Tab2GL = { data.out <- Tab2GL.conv(df,System,HZY.Red,DRB345.Flag,Cores) } )

  # Output converted file
  switch(Output,
         R = return(data.out),
         txt = write.table(data.out,file=fileName,sep="\t",quote=F,col.names=T,row.names=F),
         csv = write.csv(data.out,file=fileName,quote=F,col.names=T,row.names=F),
         pypop = write.table(data.out,file=fileName,sep="\t",quote=F,col.names=T,row.names=F) )

}



