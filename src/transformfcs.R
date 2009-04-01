#library(Biobase)
#library(prada)

transformFCS<-function(cytoFrame, filename){

 k<-readFCS(cytoFrame)
 e <- exprs(k)
 d <- description(k)

 ampliPar <- d[grep("[$]P[0-9]?E", names(d))]
 range <- as.integer(d[grep("[$]P[0-9]?R", names(d))])
 ampli <- do.call("rbind", lapply(ampliPar, function(x) as.integer(unlist(strsplit(x,",")))))
 for (i in 1:nrow(ampli)) {
   if (ampli[i, 1] > 0) {
     e[,i] <- 10^((e[,i]/(range[i] - 1)) * ampli[i, 1])  
   }
 }
 exprs(k) <- e
 return(e)
}
