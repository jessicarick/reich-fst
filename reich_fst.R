## reich fst estimator
## code by J. Rick, March 2020 -- https://github.com/jessicarick/reich-fst/
## input = genlight object
## FST will be calculated between pops in genlight object, so make sure pops are set prior to calling reich.fst()
## specify number of bootstraps using e.g. "bootstrap=100"

reich.fst <- function(gl, bootstrap=FALSE, verbose=TRUE) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    install.packages("matrixStats")
    library(matrixStats, character.only=T)
  }
  if (!require("combinat",character.only=T, quietly=T)) {
    install.packages("combinat")
  }
  if (!require("beepr",character.only=T, quietly=T)) {
    install.packages("beepr")
  }  
  
  nloc <- gl@n.loc
  npop <- length(levels(gl@pop))

  fsts <- matrix(nrow=npop,
                 ncol=npop,
                 dimnames=list(levels(gl@pop),levels(gl@pop)))
  
  if (bootstrap != FALSE){
    n.bs <- bootstrap
    bs <- data.frame(matrix(nrow=nrow(combinat::combn2(levels(gl@pop))),
                            ncol=n.bs+5))
  }
  
  k <- 0
  
  for (p1 in levels(gl@pop)){
    for (p2 in levels(gl@pop)){
      if (which(levels(gl@pop) == p1) < which(levels(gl@pop) == p2)) {
        k <- 1+k
        
        pop1 <- gl.keep.pop(gl, p1, mono.rm=FALSE, v=0)
        pop2 <- gl.keep.pop(gl, p2, mono.rm=FALSE, v=0)
        
        a1 <- colSums2(as.matrix(pop1),na.rm=T)
        a2 <- colSums2(as.matrix(pop2),na.rm=T)
        n1 <- apply(as.matrix(pop1),2,function(x) 2*sum(!is.na(x)))
        n2 <- apply(as.matrix(pop2),2,function(x) 2*sum(!is.na(x)))
        
        h1 <- (a1*(n1-a1))/(n1*(n1-1))
        h2 <- (a2*(n2-a2))/(n2*(n2-1))
        
        N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
        D <- N + h1 + h2
        
        F <- sum(N, na.rm=T)/sum(D, na.rm=T)
        fsts[p2,p1] <- F
        if (verbose == TRUE) {
        	print(paste("Pop1: ",p1,", Pop2: ",p2,", Reich FST: ",F,sep=""))
        }
        
        if (bootstrap != FALSE) {
          if (verbose == TRUE) {
            print("beginning bootstrapping")
          }
          
          bs[k,1:3] <- c(p2,p1,F)
          
          for (i in 1:n.bs){
            loci <- sample((1:nloc), nloc, replace=TRUE)
          
            pop1.bs <- matrix(as.matrix(pop1)[,loci],
                              ncol=length(loci))
            pop2.bs <- matrix(as.matrix(pop2)[,loci],
                              ncol=length(loci))
          
            a1 <- colSums2(as.matrix(pop1.bs),na.rm=T)
            a2 <- colSums2(as.matrix(pop2.bs),na.rm=T)
            n1 <- apply(as.matrix(pop1.bs),2,function(x) 2*sum(!is.na(x)))
            n2 <- apply(as.matrix(pop2.bs),2,function(x) 2*sum(!is.na(x)))
          
            h1 <- (a1*(n1-a1))/(n1*(n1-1))
            h2 <- (a2*(n2-a2))/(n2*(n2-1))
          
            N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
            D <- N + h1 + h2
          
            F.bs <- sum(N, na.rm=T)/sum(D, na.rm=T)
            bs[k,i+5] <- F.bs
          }
          if (verbose == TRUE){
            print(paste("bootstrapping 95% CI: ",
                        quantile(bs[k,6:(n.bs+5)],c(0.025),na.rm=T),"-",
                        quantile(bs[k,6:(n.bs+5)],c(0.975),na.rm=T)))
          }
          
          bs[k,4:5] <- c(quantile(bs[k,6:n.bs+5],c(0.025),na.rm=T),
                         quantile(bs[k,6:n.bs+5],c(0.975),na.rm=T))
        }
        
      }
    }
  }
  
  fsts[fsts < 0] <- 0
  colnames(bs)[1:5] <- c("pop1","pop2","fst_estimate","95_min_CI","95_max_CI")
  
  if (bootstrap != FALSE){
    fst.list <- list(fsts,bs)
    names(fst.list) <- c("fsts","bootstraps")
  } else {
    fst.list <- list(fsts)
    names(fst.list) <- "fsts"
  }
  
  return(fst.list)
  beepr::beep()
}
