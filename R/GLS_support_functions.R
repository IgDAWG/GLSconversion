#' File Name Extraction
#'
#' Function to extract file path.
#' @param x File name.
getPath <- function(x) {

  tmp.list <- list()

  if( grepl("/",x) ) {
    tmp <- unlist(strsplit(x,split="/"))
    tmp <- tmp[-which(tmp=="")]
    tmp.list$name <- tmp[length(tmp)]
    tmp.list$path <- paste("/",paste(tmp[seq(1,length(tmp)-1)],collapse="/"),sep="")
  } else {
    tmp.list$name <- x
    tmp.list$path <- NA
  }

  return(tmp.list)

}

#' Even Number Determination
#'
#' Function to identify even numbers.
#' @param x Number in question.
is.even <- function(x) {

  return(x %% 2 == 0)

}
