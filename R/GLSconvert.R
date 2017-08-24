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
GLSconvert <- function(Data,Convert,Output="txt",System="HLA",HZY.Red=FALSE,DRB345.Flag=FALSE,Cores.Lim=1L) {

  # Check Parameters
  if(is.na(match(Convert,c("GL2Tab","Tab2GL")))) { Err.Log("P.Convert") ; stop("Conversion Stopped.",call.=FALSE) }
  if(is.na(match(Output,c("R","txt","csv","pypop")))) { Err.Log("P.Output") ; stop("Conversion Stopped.",call.=FALSE) }
  if(is.na(match(System,c("HLA","KIR")))) { Err.Log("P.System") ; stop("Conversion Stopped.",call.=FALSE) }
  if(!is.logical(HZY.Red)) { Err.Log("P.HZY") ; stop("Conversion Stopped.",call.=FALSE) }
  if(!is.logical(DRB345.Flag)) { Err.Log("P.DRB") ; stop("Conversion Stopped.",call.=FALSE) }
  if(!is.numeric(Cores.Lim) || !is.integer(Cores.Lim)) { Err.Log("P.Cores") ; stop("Conversion Stopped.",call.=FALSE) }

  # MultiCore Limitations
  if (Cores.Lim!=1L) {
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
    } else { Err.Log(File.Error,Data) ; stop("Conversion Stopped.",call.=FALSE) }
  } else { df <- Data ; fileName <- "Converted.txt" }
  df[] <- lapply(df, as.character)


  # Run Conversion
  switch(Convert,
         GL2Tab = { data.out <- GL2Tab.wrapper(df,System,DRB345.Flag,Cores) } ,
         Tab2GL = { data.out <- Tab2GL.wrapper(df,System,HZY.Red,DRB345.Flag,Cores) } )

  # Output converted file
  switch(Output,
         R = return(data.out),
         txt = write.table(data.out,file=fileName,sep="\t",quote=F,col.names=T,row.names=F),
         csv = write.csv(data.out,file=fileName,quote=F,col.names=T,row.names=F),
         pypop = write.table(data.out,file=fileName,sep="\t",quote=F,col.names=T,row.names=F) )

}
