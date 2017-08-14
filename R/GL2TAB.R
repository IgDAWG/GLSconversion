GL2Tab.wrapper <- function(df,Abs.Allele) {


  lapply()


}

GL2Tab.conv <- function(x,Loci,Loci.Grp) {

  # Break GL String
  tmp <- unlist(strsplit(x,"\\^"))
  tmp <- sapply(tmp,FUN=function(x)strsplit(x,"\\+"))
  Calls <- unlist(tmp) ; names(Calls) <- NULL

  # Get Loci and Initialize Table
  Loci <- unique(unlist(lapply(strsplit(Calls,"\\*"),"[",1)))
  Loci.Grp <- rep(Loci,each=2)
  Tab <- mat.or.vec(nr=1,nc=length(Loci.Grp)) ; colnames(Tab) <- Loci.Grp
  colnames(Tab)[seq(1,length(Loci.Grp),by=2)] <- paste(Loci,"_1")
  colnames(Tab)[seq(2,length(Loci.Grp),by=2)] <- paste(Loci,"_2")

  # Populate Table
  for(i in Loci) {

    getCalls <- grep(i,Calls)

    # Assumptions for DRB345
    if(i=="HLA-DRB3" || i=="HLA-DRB4" || i=="HLA-DRB5") {

      DRB345.zygosity(i,Calls[grep("DRB",Calls)])

      if(length(getCalls)>1) {
        # For Heterzygous DRB345

        Tab[1,grep(i,colnames(Tab))] <- Calls[getCalls]
      } else {
        # For possible homozygoues DRB345



        Tab[1,grep(i,colnames(Tab))[1]] <- Calls[getCalls]
      }
    } else {
      Tab[1,grep(i,colnames(Tab))] <- Calls[grep(i,Calls)]
    }

  }# End i

}
