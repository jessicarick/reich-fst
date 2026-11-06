###########
## just the core fst calculation
calc.reich.fst <- function(gl, loc=gl@loc.names) {
  pop1 <- gl.keep.loc(gl.keep.pop(gl, levels(gl@pop)[1], mono.rm=FALSE, v=0), loc, v=0)
  pop2 <- gl.keep.loc(gl.keep.pop(gl, levels(gl@pop)[2], mono.rm=FALSE, v=0), loc, v=0)
  
  a1 <- glSum(pop1)
  a2 <- glSum(pop2)
  n1 <- -(adegenet::glNA(pop1)) + 2*nInd(pop1)
  n2 <- -(adegenet::glNA(pop2)) + 2*nInd(pop2)
  
  h1 <- (a1*(n1-a1))/(n1*(n1-1))
  h2 <- (a2*(n2-a2))/(n2*(n2-1))
  
  N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
  D <- N + h1 + h2
  
  F <- sum(N, na.rm=T)/sum(D, na.rm=T)
}

###########
## reich fst estimator from gl object
## vectorized version
## input=genlight object
## FST will be calculated between pops in genlight object
## specify number of bootstraps using "bootstrap=100"

reich.fst <- function(gl, bootstrap=FALSE, plot=FALSE, verbose=TRUE) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    exit("Function requires 'matrixStats' package; please install.")
  }
  if (!require("dplyr",character.only=T, quietly=T)) {
    exit("Function requires 'dplyr' package; please install.")
  }
  if (!require("dartR",character.only=T, quietly=T) & !require("dartRverse",character.only=T, quietly=T)) {
    exit("Function requires 'dartR' or 'dartRverse' package; please install")
  }
  if (!require("combinat",character.only=T, quietly=T)) {
    exit("Function requires 'combinat' packages; please install.")
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
        
        a1 <- glSum(pop1)
        a2 <- glSum(pop2)
        n1 <- -(adegenet::glNA(pop1)) + 2*nInd(pop1)
        n2 <- -(adegenet::glNA(pop2)) + 2*nInd(pop2)
        
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

###########
# locus-specific fst calculation
# CAUTION: WILL BE SLOW IF MANY LOCI
# genlight object needs to have only two populations

loc.reich.fst <- function(gl, pop1, pop2, 
                          loc_subset=FALSE, 
                          plot=FALSE) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    exit("Function requires 'matrixStats' package; please install.")
  }
  if (!require("dplyr",character.only=T, quietly=T)) {
    exit("Function requires 'dplyr' package; please install.")
  }
  if (!require("dartR",character.only=T, quietly=T) & !require("dartRverse",character.only=T, quietly=T)) {
    exit("Function requires 'dartR' or 'dartRverse' package; please install")
  }
  
  if(loc_subset != FALSE) {
    print("subsetting to loci specified, n=", length(loc_subset))
    gl <- gl.keep.loc(gl, loc.list=loc_subset)
  } else {
    print("no subset of loci specified, running on all loci.")
  }
  
  nloc <- gl@n.loc
  
  if(pop1 %in% levels(gl@pop) & pop2 %in% levels(gl@pop)) {
    print("both populations found in genlight object; proceeding with fst calculation.")
  } else {
    stop("one or both populations not found in genlight object! stopping now.")
  }
  # npop <- length(levels(gl@pop))
  # if(npop != 2) {
  #   stop("number of populations is not 2! this function will not work.")
  # }
  
  fsts <- tibble(loc=gl@loc.names,
                 chr=gl@chromosome,
                 pos=gl@position,
                 fst=numeric(length(gl@loc.names)))
  
  p1 <- pop1
  p2 <- pop2
  
  message(paste0("Calculating locus-specific FST values on ",length(gl@loc.names)," loci"))
  message(paste0("between populations ",p1," and ",p2))
  progress_bar = txtProgressBar(min=0, max=length(gl@loc.names), style = 1, char="*")
  
  # for(i in 1:length(gl@loc.names)){
  #   #print(paste("locus ",k," of ",length(gl@loc.names)))
  #   F <- calc.reich.fst(gl,loc=gl@loc.names[i])
  #   fsts[i,4] <- F
  #   
  #   setTxtProgressBar(progress_bar, value = i)
  # }      
  
  pop1_gl <- gl.keep.pop(gl, p1, mono.rm=FALSE, v=0)
  pop2_gl <- gl.keep.pop(gl, p2, mono.rm=FALSE, v=0)
  
  a1 <- glSum(pop1_gl)
  a2 <- glSum(pop2_gl)
  n1 <- -(adegenet::glNA(pop1_gl)) + 2*nInd(pop1_gl)
  n2 <- -(adegenet::glNA(pop2_gl)) + 2*nInd(pop2_gl)
  
  h1 <- (a1*(n1-a1))/(n1*(n1-1))
  h2 <- (a2*(n2-a2))/(n2*(n2-1))
  
  N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
  D <- N + h1 + h2
  
  F <- N / D
  
  # make negative values 0
  F[F < 0] <- 0
  
  # add to tibble
  fsts$fst <- F
  
  if (plot == TRUE){
    message("drawing plot of locus-specific FST estimates")
    
    if (!require("ggplot2",character.only=T, quietly=T)) {
      install.packages("ggplot2")
      library(ggplot2, character.only=T)
    }
    
    fst.plot <- ggplot(fsts, aes(x=loc,y=fst)) + 
      geom_point(size=1, alpha=0.6) +
      geom_hline(yintercept=mean(fsts$fst,na.rm=T), lty=3, lwd=0.5, col="gray50") +
      theme_bw() +
      xlab("Locus") +
      ylab("Reich-Patterson FST Estimate") +
      theme(legend.position="none",
            axis.text.x=element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank())
    
    print(fst.plot)
  }
  
  return(fsts)
}

###########
# windowed fst calculation
# CAUTION: WILL BE SLOW IF MANY LOCI
# genlight object needs to have only two populations

win.reich.fst <- function(gl, plot=FALSE, chrom=TRUE,
                          win.type="snp", win.size=100) { 
  if (!require("matrixStats",character.only=T, quietly=T)) {
    exit("Function requires 'matrixStats' package; please install.")
  }
  if (!require("dplyr",character.only=T, quietly=T)) {
    exit("Function requires 'dplyr' package; please install.")
  }
  if (!require("dartR",character.only=T, quietly=T) & !require("dartRverse",character.only=T, quietly=T)) {
    exit("Function requires 'dartR' or 'dartRverse' package; please install")
  }
  
  nloc <- gl@n.loc
  npop <- length(levels(gl@pop))
  if(npop != 2) {
    stop("number of populations is not 2! this function will not work.")
  }
  
  message(paste0("Calculating windowed FST values on ",length(gl@loc.names)," loci"))
  message(paste0("between populations ",levels(gl@pop)[1]," and ",levels(gl@pop)[2]))
  
  # creating windows
  if (win.type == "snp") {
    if (chrom == TRUE) {
      message(paste("working with",win.size,"SNP windows and separating by chromosome"))
      
      ## first calculate fst for each snp
      pop1_gl <- gl.keep.pop(gl, pop1, mono.rm=FALSE, v=0)
      pop2_gl <- gl.keep.pop(gl, pop2, mono.rm=FALSE, v=0)
      
      a1 <- glSum(pop1_gl)
      a2 <- glSum(pop2_gl)
      n1 <- -(adegenet::glNA(pop1_gl)) + 2*nInd(pop1_gl)
      n2 <- -(adegenet::glNA(pop2_gl)) + 2*nInd(pop2_gl)
      
      h1 <- (a1*(n1-a1))/(n1*(n1-1))
      h2 <- (a2*(n2-a2))/(n2*(n2-1))
      
      N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
      D <- N + h1 + h2
      
      ## now aggregate in windows
      windows <- tibble(win_num = integer(),
                        chrom = character(),
                        bin_start = integer(),
                        bin_end = integer(),
                        bin_mid = integer(),
                        reich_fst = numeric(),
                        n_snps = integer())
      snps <- tibble(chrom = gl@chromosome,
                     pos = gl@position,
                     loc_name = gl@loc.names,
                     N = N,
                     D = D) %>%
        group_by(chrom) %>%
        mutate(win_num = ceiling(row_number()/win.size))
      win_fsts <- snps %>%
        group_by(chrom, win_num) %>%
        summarize(bin_start = min(pos),
                  bin_end = max(pos),
                  n_snp = n(),
                  N_sum = sum(N, na.rm=T),
                  D_sum = sum(D, na.rm=T)) %>%
        mutate(bin_mid = round((bin_start + bin_end)/2,0),
               reich_fst = N_sum / D_sum) %>%
        select(chrom, win_num, bin_start, bin_end, bin_mid, reich_fst, n_snp) %>%
        group_by(chrom) %>%
        arrange(bin_start, .by_group = TRUE) %>%
        ungroup() %>%
        mutate(win_seq = seq(1,n()))
      
      ## set negative numbers to 0
      win_fsts$reich_fst[win_fsts$reich_fst < 0] <- 0
      
      # if (min(n_win$n_win) < 5) {
      #   message("WARNING: Some chromosomes have fewer than 5 windows; consider using a smaller window size.")
      # } else {
      #   message("All chromosomes have at least 5 windows.")
      # }
      
      # k <- 1
      # progress_bar = txtProgressBar(min=0, max=sum(n_win$n_win), style = 1, char="*")
      # 
      # for (chr in unique(snps$chrom)){
      #   snps_chr <- snps %>% filter(chrom == chr)
      #   
      #   for (win in unique(snps_chr$win_num)){
      #     fst_win <- calc.reich.fst(gl, loc=snps_chr %>% filter(win_num == win) %>% pull(loc_name))
      #     
      #     windows <- windows %>%
      #       add_row(win_num = k,
      #               chrom = chr,
      #               bin_start = min(snps_chr %>% filter(win_num == win) %>% pull(pos)),
      #               bin_end = max(snps_chr %>% filter(win_num == win) %>% pull(pos)),
      #               bin_mid = floor(bin_start + (bin_end-bin_start)/2),
      #               reich_fst = fst_win,
      #               n_snps = nrow(snps_chr %>% filter(win_num == win)))
      #     
      #     k <- k+1
      #     setTxtProgressBar(progress_bar, value = k)
      #   }
      #   
      # }
      
    } else {
      message(paste("working with",win.size,"SNP windows and NOT separating by chromosome"))
      
      ## first calculate fst for each snp
      pop1_gl <- gl.keep.pop(gl, pop1, mono.rm=FALSE, v=0)
      pop2_gl <- gl.keep.pop(gl, pop2, mono.rm=FALSE, v=0)
      
      a1 <- glSum(pop1_gl)
      a2 <- glSum(pop2_gl)
      n1 <- -(adegenet::glNA(pop1_gl)) + 2*nInd(pop1_gl)
      n2 <- -(adegenet::glNA(pop2_gl)) + 2*nInd(pop2_gl)
      
      h1 <- (a1*(n1-a1))/(n1*(n1-1))
      h2 <- (a2*(n2-a2))/(n2*(n2-1))
      
      N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
      D <- N + h1 + h2
      
      ## now aggregate in windows
      snps <- tibble(chrom = gl@chromosome,
                     pos = gl@position,
                     loc_name = gl@loc.names,
                     N = N,
                     D = D) %>%
        #group_by(chrom) %>%
        mutate(win_num = ceiling(row_number()/win.size))
      win_fsts <- snps %>%
        group_by(win_num) %>%
        summarize(bin_start = min(pos),
                  bin_end = max(pos),
                  n_snp = n(),
                  N_sum = sum(N, na.rm=T),
                  D_sum = sum(D, na.rm=T)) %>%
        mutate(bin_mid = round((bin_start + bin_end)/2,0),
               reich_fst = N_sum / D_sum) %>%
        select(win_num, bin_start, bin_end, bin_mid, reich_fst, n_snp)
      
      ## set negative numbers to 0
      win_fsts$reich_fst[win_fsts$reich_fst < 0] <- 0
      
    }
    
  } else if (win.type == "bp") {
    if (chrom == TRUE) {
      message(paste("working with",win.size,"bp windows and separating by chromosome"))
    
      ## first calculate fst for each snp
      pop1_gl <- gl.keep.pop(gl, pop1, mono.rm=FALSE, v=0)
      pop2_gl <- gl.keep.pop(gl, pop2, mono.rm=FALSE, v=0)
      
      a1 <- glSum(pop1_gl)
      a2 <- glSum(pop2_gl)
      n1 <- -(adegenet::glNA(pop1_gl)) + 2*nInd(pop1_gl)
      n2 <- -(adegenet::glNA(pop2_gl)) + 2*nInd(pop2_gl)
      
      h1 <- (a1*(n1-a1))/(n1*(n1-1))
      h2 <- (a2*(n2-a2))/(n2*(n2-1))
      
      N <- (a1/n1 - a2/n2)^2 - h1/n1 - h2/n2
      D <- N + h1 + h2
      
      ## now aggregate in windows
      snps <- tibble(chrom = gl@chromosome,
                     pos = gl@position,
                     loc_name = gl@loc.names,
                     N = N,
                     D = D) %>%
        group_by(chrom) %>%
        mutate(win_num = ceiling(pos/win.size))
      win_fsts <- snps %>%
        group_by(chrom, win_num) %>%
        summarize(bin_start = min(pos),
                  bin_end = max(pos),
                  n_snp = n(),
                  N_sum = sum(N, na.rm=T),
                  D_sum = sum(D, na.rm=T)) %>%
        mutate(bin_mid = round((bin_start + bin_end)/2,0),
               reich_fst = N_sum / D_sum) %>%
        select(chrom, win_num, bin_start, bin_end, bin_mid, reich_fst, n_snp) %>%
        group_by(chrom) %>%
        arrange(bin_start, .by_group = TRUE) %>%
        ungroup() %>%
        mutate(win_seq = seq(1,n()))
      
      ## set negative numbers to 0
      win_fsts$reich_fst[win_fsts$reich_fst < 0] <- 0
      
      
    } else {
      stop("it doesn't make sense to work with bp positions and not separate by chromosome! try again.")
      
    }
    
  } else {
    stop("unrecognized win.type argument, should be either snp or bp")
  }
  
  # make negative values 0
  windows$reich_fst[windows$reich_fst < 0] <- 0
  
  if (plot == TRUE & chrom == TRUE){
    message("drawing plot of windowed FST estimates")
    
    if (!require("ggplot2",character.only=T, quietly=T)) {
      install.packages("ggplot2")
      library(ggplot2, character.only=T)
    }
    
    fst.plot <- ggplot(win_fsts, aes(x=win_seq,y=reich_fst)) + 
      geom_point(aes(color=chrom), size=0.5, alpha=0.6) +
      geom_hline(yintercept=mean(win_fsts$reich_fst,na.rm=T), lty=3, lwd=0.5, col="gray50") +
      theme_bw() +
      scale_color_manual(values=rep(c("gray60","gray20"),nrow(win_fsts))) +
      xlab("Locus") +
      ylab("Windowed Reich-Patterson FST") +
      theme(legend.position="none",
            axis.text.x=element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank())
    
    print(fst.plot)
  } else if (plot == TRUE & chrom == FALSE) {
    message("drawing plot of windowed FST estimates")
    
    if (!require("ggplot2",character.only=T, quietly=T)) {
      install.packages("ggplot2")
      library(ggplot2, character.only=T)
    }
    
    fst.plot <- ggplot(win_fsts, aes(x=win_seq,y=reich_fst)) + 
      geom_point(aes(color=chrom), size=0.5, alpha=0.6) +
      geom_hline(yintercept=mean(win_fsts$reich_fst,na.rm=T), lty=3, lwd=0.5, col="gray50") +
      theme_bw() +
      scale_color_manual(values=rep(c("gray60","gray20"),nrow(win_fsts))) +
      xlab("Locus") +
      ylab("Windowed Reich-Patterson FST") +
      theme(legend.position="none",
            axis.text.x=element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank())
    
    print(fst.plot)
  }
  
  return(win_fsts)
}
