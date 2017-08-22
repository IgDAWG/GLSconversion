#' File Name Extraction
#'
#' Function to extract file path.
#' @param x File name.
getName <- function(x) {

  tmpDir <- dirname(x)
  tmpName <- basename(x)

  if(basename(x)==x) {
    outName <- paste("Converted_",x,sep="")
  } else {
    outName <- paste(tmpDir,"/Converted_",tmpName,sep="")
  }

  return(outName)

}
