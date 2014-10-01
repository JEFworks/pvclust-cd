#pvclust-cd
##Hierarchical Clustering with P-Values via Multiscale Bootstrap using Custom Distance Functions

## Description
This is a modification of Ryota Suzuki and Hidetoshi Shimodaira's Hierarchical Clustering with P-Values via Multiscale Bootstrap Resampling package pvclust to accept custom distance functions as opposed to being restricted to "correlation", "abscor", and "uncentered". For more information on the pvclust package, see: http://www.is.titech.ac.jp/~shimo/prog/pvclust/

## Tutorial

### Loading

    source('pvclust-cd.R')

### Creating custom distance functions

    ## note custom distance functions must take 2 vectors (not matrices)
    frac.shared <- function(x,y) {
        pseudo <- 10e-6
	maxvar <- max(sum(x), sum(y))
        fs <- 1-sum(x*y)/(maxvar+pseudo)
        if(fs>1) { fs <- 1 }
        fs
    }

    ## simulate mutations for 50 cells
    ## simulate some germline mutations
    germline.muts <- matrix(1, ncol=50, nrow=10)
    ## simulate some subpopulation specific mutations (subclonal somatic mutations)
    subpop.muts <- rbind(  c(rep(1, 10),rep(0, 40)),
            c(rep(1, 20),rep(0, 30)),
            c(rep(1, 30),rep(0, 20)),
            c(rep(1, 40),rep(0, 10))
    )
    ## simulate random mutations
    random.muts <- matrix(rbinom(100*50, 1, 0.1), ncol = 50, nrow = 50)

    ## cluster cells on just the subpopulation specific mutations 
    mat <- subpop.muts
    colnames(mat) <- paste('cell', 1:ncol(mat))
    rownames(mat) <- paste('mut', 1:nrow(mat))
    hc1 <- cd.pvclust(mat, method.hclust="ward", method.dist=frac.shared, nboot=100)
    plot(hc1)
    ## subclonal architecture should be recapitulated confidently

    ## cluster cells on all mutations
    ## ie. add noise
    mat <- rbind(germline.muts, subpop.muts, random.muts)
    colnames(mat) <- paste('cell', 1:ncol(mat))
    rownames(mat) <- paste('mut', 1:nrow(mat))
    hc2 <- cd.pvclust(mat, method.hclust="ward", method.dist=frac.shared, nboot=100)
    plot(hc2)
    ## how does noise impact the inferred subclonal architecture?
    ## is there a better custom distance function that would be more robust to noise?

### Compared to pvclust
    
    ## compare with original pvclust
    hc1 <- pvclust(t(mtcars), method.hclust="complete", method.dist="correlation", use.cor="complete.obs", nboot=10)
    plot(hc1)

    ## check that results are same
    cor.function <- function(x,y) {
        1-cor(x,y, method="pearson", use="complete.obs")
    }
    hc2 <- cd.pvclust(t(mtcars), method.hclust="complete", method.dist=cor.function, nboot=10)
    plot(hc2)
    ## much slower because correlations are computed one by one instead of as a matrix
    ## results differ slightly due to stochastic nature of calculations
