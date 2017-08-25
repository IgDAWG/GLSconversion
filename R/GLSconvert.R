#' Genotype List String Conversion
#'
#' Main Workhorse wrapper for cross converting columnar table to GL string representaion.
#' @param Data String File name or R Data Frame.
#' @param Convert String Direction for conversion.
#' @param Output String Type of output.
#' @param System String Genetic system (HLA or KIR) of the data being converted
#' @param HZY.Red Logical Reduction of homozygote genotypes to single allele.
#' @param DRB345.Check Logical Check DR haplotypes for consistency and flag unusual haplotypes.
#' @param Cores.Lim Integer How many cores can be used.
GLSconvert <- function(Data,Convert,Output="txt",System="HLA",HZY.Red=FALSE,DRB345.Check=TRUE,Cores.Lim=1L) {

  # Check Parameters
  Check.Params(Convert,Ouput,System,HZY.Red,DRB345.Check,Cores.Lim)

  # MultiCore Limitations
  if ( Cores.Lim!=1L ) {
    Cores.Max <- as.integer( floor( parallel::detectCores() * 0.9) )
    if(Sys.info()['sysname']=="Windows" && as.numeric(Cores.Lim)>1) {
      Err.Log("Windows.Cores") ; stop("Conversion stopped.",call. = F)
    } else if( Cores.Lim > Cores.Max ) { Cores <- Cores.Max
    } else { Cores <- Cores.Lim }
  } else { Cores <- Cores.Lim }

  # Nomenclature system
  if( System == "HLA" ) { System <- "HLA-" }

  # Read in Data and Set Output File Name
  if( is.character(Data) ) {
    if( file.exists(Data) ) {
      df <- read.table(file=Data,header=T,sep="\t",stringsAsFactors=FALSE)
      fileName <- getName(Data)
    } else { Err.Log("File.Error",Data) ; stop("Conversion Stopped.",call.=FALSE) }
  } else { df <- Data ; fileName <- "Converted.txt" }
  df[] <- lapply(df, as.character)

  # Run Conversion
  switch(Convert,
         GL2Tab = { data.out <- GL2Tab.wrapper(df,System,Cores) } ,
         Tab2GL = { data.out <- Tab2GL.wrapper(df,System,HZY.Red,Cores) } )

  # Output DRB.HapFlag for HLA data
  if( System=="HLA-" && !DRB345.Check ) { data.out <- data.out[,-grep('DRB.HapFlag',colnames(data.out))]  }

  # Output converted file
  switch(Output,
         R = return(data.out),
         txt = write.table(data.out,file=fileName,sep="\t",quote=F,col.names=T,row.names=F),
         csv = write.csv(data.out,file=fileName,quote=F,col.names=T,row.names=F),
         pypop = write.table(data.out,file=fileName,sep="\t",quote=F,col.names=T,row.names=F) )

}
