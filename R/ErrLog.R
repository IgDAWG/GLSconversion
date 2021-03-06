#' Error Code Display and Logging
#'
#' Displays error codes attributable to data formatting. Writes to log file.
#' @param x Log Code.
#' @param y Misc information relevant to error.
#' @param z Misc information relevant to error.
#' @note This function is for internal use only.
Err.Log.GLS <- function (x, y=NULL, z=NULL) {

  cat("*****ERROR!******\n")
  switch(x,

         #Parameters
         P.Missing = { Error <- paste("\nNo ",y," specified. This parameter is not optional. Please see vignette.",sep="") },
         P.Error = { Error <- paste("\nInvalid ",y," parameter. Please see vignette.",sep="") },
         Windows.Cores = { Error <- "\nYou have exceed the maximum allowable cores for Windows. Please see vignette." },

         #Notifications
         Tab.Format = { Error <- "\nThe conversion tool encountered GL string delimiters. This isn't valid data for Tab2GL converion. Please see vignette." },
         GL.Format = { Error <- "\nYour GL strings may not be properly formatted. The conversion tool did not encounter the classic GL string delimiters. Please see vignette." },
         File.Error = { Error <- paste("\nThe conversion tool could not locate a file labeled ",y," in the specificied working directory.",sep="") },
         GTYPE.Amb = { Error <- paste("\nThis appears to contain genotype list piping ('|') for genotype ambiguity strings (data rows: ",y,"). This is not supported in GLSconversion.",sep="") },
         Table.Col = { Error <- "\nThe table for Tab2GL conversion is not properly formatted, too few columns. Please see vignette." },
         Table.Pairs = { Error <- "\nThe table for Tab2GL conversion is not properly formatted, no locus column pairs encountered. Please see vignette." },
         Table.Amb = { Error <- "\nYour data has duplicate identifying information rows, perhaps due to data genotype ambiguity." },
         Locus.MultiField = { Error <- paste("\nYour GL string may be invalid. A locus cannot appear in multiple gene fields! ",z,ifelse(grepl(",",z)," appear"," appears")," in multiple fields of the GL string: ", y, ". Please see vignette.",sep="") },
         Allele.Amb.Format = { Error <- paste("\nYour GL string may be invalid. The ambiguous allele ",y," is not properly formatted. Please see vignette.", sep="") }
  )

  cat(Error,"\n")

}
