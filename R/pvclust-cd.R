######################################
# pvclust-cd.R
# @author: Jean Fan
# @date: July 23, 2014
# @description: Extension of pvclust R package to accept custom distance functions
# @usage: See example functions
######################################

source('pvclust.R')
source('pvclust-internal.R')

cd.pvclust <- function(data, method.hclust="ward", method.dist=method.dist, nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE)
{
    # data: (n,p) matrix, n-samples, p-variables
    n <- nrow(data); p <- ncol(data)

    # hclust for original data
    METHODS <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
    method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
    distance <- cd.dist.pvclust(data, method=method.dist)
    data.hclust <- hclust(distance, method=method.hclust)
    
    # multiscale bootstrap
    size <- floor(n*r)
    rl <- length(size)
    
    if(rl == 1) {
      if(r != 1.0)
        warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
      
      r <- list(1.0)
    }
    else
      r <- as.list(size/n)
    
    mboot <- lapply(r, cd.boot.hclust, data=data, object.hclust=data.hclust, nboot=nboot, method.dist=method.dist, method.hclust=method.hclust, store=store)
    
    result <- pvclust.merge(data=data, object.hclust=data.hclust, mboot=mboot)
    
    return(result)
}


cd.dist.pvclust <- function(mat, method)
{
    n <- ncol(mat)
    dist.mat <- matrix(0, ncol = n, nrow = n)
    colnames(dist.mat) <- rownames(dist.mat) <- colnames(mat)
    for(i in 1:n) {
        for(j in 1:n) {
            dist.mat[i,j] <- method(mat[,i],mat[,j])
            ## should introduce use pairwise complete, na.rm, and other options
        }
    }
    res <- as.dist(dist.mat)
    attr(res,"method") <- "custom"
    
    return(res)
}


cd.boot.hclust <- function(r, data, object.hclust, method.dist, method.hclust, nboot, store)
{ 
  n     <- nrow(data)
  size  <- round(n*r, digits=0)
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n

  pattern   <- hc2split(object.hclust)$pattern
  edges.cnt <- table(factor(pattern)) - table(factor(pattern))
  st <- list()
  
  # bootstrap start
  rp <- as.character(round(r,digits=2)); if(r == 1) rp <- paste(rp,".0",sep="")
  cat(paste("Bootstrap (r = ", rp, ")... ", sep=""))
  w0 <- rep(1,n) # equal weight
  na.flag <- 0
  
  for(i in 1:nboot){
    #if(weight && r>10) {  ## <- this part should be improved
    #  w1 <- as.vector(rmultinom(1,size,w0)) # resampled weight
    #  suppressWarnings(distance <- cd.distw.pvclust(data,w1,method=method.dist))
    #} else {
    #  smpl <- sample(1:n, size, replace=TRUE)
    #  suppressWarnings(distance  <- cd.dist.pvclust(data[smpl,],method=method.dist))
    #}
    # remove weighting option
    smpl <- sample(1:n, size, replace=TRUE)
    suppressWarnings(distance  <- cd.dist.pvclust(data[smpl,],method=method.dist))

    if(all(is.finite(distance))) { # check if distance is valid
      x.hclust  <- hclust(distance,method=method.hclust)
      pattern.i <- hc2split(x.hclust)$pattern # split
      edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
    } else {
      x.hclust <- NULL
	  na.flag <- 1
    }

    if(store)
      st[[i]] <- x.hclust
  }
  cat("Done.\n")
  # bootstrap done
  
  if(na.flag == 1)
	warning(paste("inappropriate distance matrices are omitted in computation: r = ", r), call.=FALSE)

  boot <- list(edges.cnt=edges.cnt, method.dist=method.dist, method.hclust=method.hclust, nboot=nboot, size=size, r=r, store=st)
  class(boot) <- "boot.hclust"
  
  return(boot)
}


cd.example <- function() {
    # simulate mutations for 100 cells
    # simulate some germline mutations
    germline.muts <- matrix(1, ncol=100, nrow=10)
    # simulate some somatic subpop muts
    subpop.muts <- rbind(  c(rep(1, 10),rep(0, 90)),
            c(rep(1, 20),rep(0, 80)),
            c(rep(1, 30),rep(0, 70)),
            c(rep(1, 40),rep(0, 60)),
            c(rep(1, 50),rep(0, 50)),
            c(rep(1, 60),rep(0, 40)),
            c(rep(1, 70),rep(0, 30)),
            c(rep(1, 80),rep(0, 20)),
            c(rep(1, 90),rep(0, 10)),
            c(rep(0, 70),rep(1, 30)),
            c(rep(0, 80),rep(1, 20)),
            c(rep(0, 90),rep(0, 10))
    )
    
    random.muts <- matrix(rbinom(100*50, 1, 0.1), ncol = 100, nrow = 50)
    
    frac.shared <- function(x,y) {
        maxvar <- max(sum(x), sum(y))
        1-sum(x*y)/maxvar
    }

    #frac.shared <- function(x,y) {
    #    1-sum(x*y)/length(x)
    #}
    
    mat <- rbind(germline.muts, subpop.muts)
    colnames(mat) <- paste('cell', 1:ncol(mat))
    rownames(mat) <- paste('mut', 1:nrow(mat))
    
    hc <- cd.pvclust(mat, method.hclust="ward", method.dist=frac.shared, nboot=100)
    plot(hc)
    
    
    mat <- rbind(germline.muts, subpop.muts, random.muts)
    colnames(mat) <- paste('cell', 1:ncol(mat))
    rownames(mat) <- paste('mut', 1:nrow(mat))
    
    hc <- cd.pvclust(mat, method.hclust="ward", method.dist=frac.shared, nboot=100)
    plot(hc)
    cutree(hc$hclust, 10)
    
}


cd.example2 <- function() {
    
    b62 <- as.matrix(read.table("ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62", check.names=FALSE))
    ## Blosum function to quantify the difference between two sequences
    blosum <- function(a,b) {
        score <- sum(abs(b62[a[!(a %in% b)], b[!(b %in% a)]]))
    }
    aa.seqs <- c(a="ANQGH",b="ANCGH",c="ANQEH",d="ANQES",e="RDCGH",f="RNCGH")
    mat <- t(do.call(rbind,lapply(aa.seqs, function(s) strsplit(s, "")[[1]])))
    
    hc <- cd.pvclust(mat, method.hclust="complete", method.dist=blosum, nboot=10)
    plot(hc)

}

cd.example3 <- function() {
    
    ## check that results are same
    cor.function <- function(x,y) {
        1-cor(x,y, method="pearson", use="complete.obs")
    }
    hc <- cd.pvclust(t(mtcars), method.hclust="complete", method.dist=cor.function, nboot=10)
    plot(hc)
    ## results are the same but much slower because correlations are computed one by one instead of as a matrix
    
    ## using built in
    hc2 <- pvclust(t(mtcars), method.hclust="complete", method.dist="correlation", use.cor="complete.obs", nboot=10)
    plot(hc2)

}