## reich fst estimator
## vectorized version
## input=genlight object
## FST will be calculated between pops in genlight object
## specify number of bootstraps using "bootstrap=100"

reich.fst <- function(gl, bootstrap=FALSE, plot=FALSE, verbose=TRUE) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    install.packages("matrixStats")
    library(matrixStats, character.only=T)
  }
  if (!require("dplyr",character.only=T, quietly=T)) {
    install.packages("dplyr")
    library(dplyr, character.only=T)
  }
  if (!require("dartR",character.only=T, quietly=T)) {
    install.packages("dartR")
    library(dartR, character.only=T)
  }
  if (!require("combinat",character.only=T, quietly=T)) {
    install.packages("combinat")
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
          
          bs[k,1:3] <- c(p2,p1,as.numeric(F))
          
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
                        quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),"-",
                        quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T)))
          }
          
          bs[k,4:5] <- c(quantile(bs[k,6:(n.bs+5)],0.025,na.rm=T),
                         quantile(bs[k,6:(n.bs+5)],0.975,na.rm=T))
        }
        
      }
    }
  }
  
  fsts[fsts < 0] <- 0
  
  if (bootstrap != FALSE){
    colnames(bs)[1:5] <- c("pop1","pop2","fst_estimate","min_CI","max_CI")
    fst.list <- list(fsts,bs)
    names(fst.list) <- c("fsts","bootstraps")
    
    if (plot == TRUE){
      print("drawing plot with bootstraps")
      
      if (!require("ggplot2",character.only=T, quietly=T)) {
        install.packages("ggplot2")
        library(ggplot2, character.only=T)
      }
      
      plot.data <- bs[,1:5]
      plot.data$fst_estimate <- as.numeric(plot.data$fst_estimate)
      plot.data$min_CI <- as.numeric(plot.data$min_CI)
      plot.data$max_CI <- as.numeric(plot.data$max_CI)
      plot.data$pop_pair <- paste(plot.data$pop1,plot.data$pop2,sep="_")
      plot.data$signif <- case_when(plot.data$min_CI > 0 ~ TRUE,
                                    TRUE ~ FALSE)

      
      bs.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate,col=signif)) + 
        geom_point(size=2) + 
        coord_flip() + 
        geom_errorbar(aes(ymin=min_CI,ymax=max_CI),width=0.1,linewidth=1) + 
        geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
        theme_minimal() + 
        theme(legend.position="none")
      
      print(bs.plot)
    }
  } else {
    fst.list <- list(fsts)
    names(fst.list) <- "fsts"
    
    if (plot == TRUE){
      print("drawing plot without bootstraps")
      
      if (!require("ggplot2",character.only=T, quietly=T)) {
        install.packages("ggplot2")
        library(ggplot2, character.only=T)
      }
      
      plot.data <- data.frame(combinat::combn2(row.names(fsts)),
                              fst_estimate=fsts[lower.tri(fsts)])
      plot.data$pop_pair <- paste(plot.data$X1,plot.data$X2,sep="_")

      fst.plot <- ggplot(plot.data, aes(x=pop_pair,y=fst_estimate)) + 
        geom_point(size=2) + 
        coord_flip() + 
        geom_hline(yintercept=0, lty=2, lwd=1, col="gray50") + 
        theme_minimal() + 
        theme(legend.position="none")
      
      print(fst.plot)
    }
  }
  
  return(fst.list)
  beepr::beep()
}

# locus-specific fst calculation
# CAUTION: WILL BE SLOW IF MANY LOCI
# genlight object needs to have only two populations

loc.reich.fst <- function(gl, plot=FALSE, verbose=TRUE) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    install.packages("matrixStats")
    library(matrixStats, character.only=T)
  }
  if (!require("dplyr",character.only=T, quietly=T)) {
    install.packages("dplyr")
    library(dplyr, character.only=T)
  }
  
  nloc <- gl@n.loc
  npop <- length(levels(gl@pop))
  if(npop != 2) {
    stop("number of populations is not 2! this function will not work.")
  }
  
  fsts <- tibble(loc=gl@loc.names,
                 chr=gl@chromosome,
                 pos=gl@position,
                 fst=numeric(length(gl@loc.names)))
  
  p1 <- levels(gl@pop)[1]
  p2 <- levels(gl@pop)[2]
  
  print(paste0("Calculating locus-specific FST values on ",length(gl@loc.names)," loci"))
  print(paste0("between populations ",p1," and ",p2))
  progress_bar = txtProgressBar(min=0, max=length(gl@loc.names), style = 1, char="*")
  
  for(i in 1:length(gl@loc.names)){
        #print(paste("locus ",k," of ",length(gl@loc.names)))
        pop1 <- gl.keep.pop(gl.keep.loc(gl,gl@loc.names[i],v=0), p1, mono.rm=FALSE, v=0)
        pop2 <- gl.keep.pop(gl.keep.loc(gl,gl@loc.names[i],v=0), p2, mono.rm=FALSE, v=0)
        
        a1 <- colSums2(as.matrix(pop1),na.rm=T)
        a2 <- colSums2(as.matrix(pop2),na.rm=T)
        n1 <- apply(as.matrix(pop1),2,function(x) 2*sum(!is.na(x)))
        n2 <- apply(as.matrix(pop2),2,function(x) 2*sum(!is.na(x)))
        
        h1 <- (a1*(n1-a1))/(n1*(n1-1))
        h2 <- (a2*(n2-a2))/(n2*(n2-1))
        
        N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
        D <- N + h1 + h2
        
        F <- sum(N, na.rm=T)/sum(D, na.rm=T)
        fsts[i,4] <- F
        
        setTxtProgressBar(progress_bar, value = i)
  }      
  
  # make negative values 0
  fsts$fst[fsts$fst < 0] <- 0
  
  return(fsts)
}
