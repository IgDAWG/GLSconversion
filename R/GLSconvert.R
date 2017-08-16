#' Genotype List String Conversion
#'
#' Main Workhorse wrapper for cross converting columnar table to GL string representaion.
#' @param Data String File name or data frame.
#' @param Conv String Direction for conversion.
#' @param Output String Type of output.
#' @param HZY.Red Logical Reduction of homozygote genotypes to single allele.
#' @param DRB345.Check Logical Flag unusual DR haplotypes.
#' @param Cores.Lim Integer How many cores can be used.
GLS.Convert <- function(Data,Conv="GL2Tab",Output="pypop",HZY.Red=TRUE,DRB345.Check=TRUE,Cores.Lim=1L) {

  # All input files must be tab delimited text files of unamibugous HLA data.
  # All genotype calls should be formatted as loci*allele format (e.g., HLA-A*01:01:01:01)
  # This processing script makes certain assumptions based on Class I and Class II loci as well as
  #   known DR haplotypes.
  # Data:
  #    R data frame object or tab delimited text file name (full path recommended)
  # Conv: Which direction for conversion?
  #    GL2Tab - Expand GL string data (3 column) to tabular format (default)
  #    Tab2GL - Condense tabular data to GL string format
  # Output: what is output type
  #    R - R data frame object
  #    txt - outputs to Tab delimited text file (default)
  #    csv - outputs to CSV delimited text file
  #    pypop - same as txt with *.pop file extension
  # ....................................................... Paramters with HLA class I and class II assumptions
  # HZY.red:
  #    Homozygous reduction: Should non-dRB345 homozygotes be represent by a single allele name in GL string?
  #    For example: HLA-A*01:01:01:01 + HLA-A*01:01:01:01 as HLA-A*01:01:01:01
  #    Default behavior is to reduce HLA genotypes for homozygotes in the GL string
  #    This setting does not impact DRB3, DRB4, or DRB5 genotype calls
  # DRB345.Check
  #    Should DR haplotypes be assessed for correct zygosity and unusual DR haplotypes be flagged with a '!'
  #    appended to genotype?

  library(parallel)

  # MultiCore Limitations
  if (Cores.Lim!=1L) {
    Cores.Max <- as.integer( floor( parallel::detectCores() * 0.9) )
    if(Sys.info()['sysname']=="Windows" && as.numeric(Cores.Lim)>1) {
      stop("Windows cores max exceeded. Conversion stopped.",call. = F)
    } else if( Cores.Lim > Cores.Max ) { Cores <- Cores.Max
    } else { Cores <- Cores.Lim }
  } else { Cores <- Cores.Lim }


  if(is.character(Data)) {
    df <- read.table(file=Data,header=T,sep="\t",stringsAsFactors=FALSE)
  } else { df <- Data }

  switch(Conv,
         GL2Tab = { data.out <- GL2Tab.conv(df,DRB345.Check,Cores.Lim) } ,
         Tab2GL = { data.out <- TAB2GL.conv(df,HZY.Red,Output,Cores.Lim) } )

  switch(Output,
         R = return(data.out),
         txt = write.table(data.out,file="Converted.txt",sep="\t",quote=F,col.names=T,row.names=F),
         csv = write.csv(data.out,file="Converted.csv",quote=F,col.names=T,row.names=F),
         pypop = write.table(data.out,file="Converted.pop",quote=F,col.names=T,row.names=F) )

}



