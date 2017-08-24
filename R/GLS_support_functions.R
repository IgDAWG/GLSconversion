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


#' Error Code Display and Logging
#'
#' Displays error codes attributable to data formatting. Writes to log file.
#' @param x Log Code.
#' @param y Misc information relevant to error.
#' @note This function is for internal BIGDAWG use only.
Err.Log <- function (x, y=NULL) {

  switch(x,
         #Parameters
         P.Convert = { Error <- "\nInvalid Convert parameter. Please see vignette." },
         P.Output =  { Error <- "\nInvalid Output parameter. Please see vignette." },
         P.System =  { Error <- "\nInvalid System parameter. Please see vignette." },
         P.HZY = { Error <- "\nInvalid HZY.Red parameter. Please see vignette." },
         P.DRB = { Error <- "\nInvalid DRB345.Flag parameter. Please see vignette." },
         P.Cores = { Error <- "\nInvalid Cores.Lim parameter. Please see vignette." },
         Windows.Cores = { Error <- "\nYou have exceed the maximum allowable cores for Windows. Please see vignette." },
         #Notifications
         File.Error = { Error <- paste("\nThe conversion tool could not locate a file labeled ",y," in the specificied working directory.",sep="") },
         GTYPE.Amb = { Error <- "\nThis appears to contain genotype list piping ('|') for genotype ambiguity strings. This is not supported in GLSconversion." },
         Table.Col = { Error <- "\nThe table for Tab2GL conversion is not properly formatted, too few columns. Please see vignette." },
         Table.Amb = { Error <- "\nYour data has duplicate identifying information rows, perhaps due to data genotype ambiguity." }
  )

  cat(Error,"\n")

}
