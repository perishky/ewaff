## create variable y, such that cor(x,y) == r
rcor <- function(x,r) {
    e <- rnorm(length(x), mean=0, sd=sqrt(1-r^2))
    mx <- mean(x)
    sx <- sd(x)
    s <- (x-mx)/sx
    (r*s + e)*sx + mx
}

## clusterMaker
## 
## Code in this file from:
## 
## Jaffe AE, Murakami P, Lee H, Leek JT, Fallin DM, Feinberg AP, Irizarry
## RA (2012). Bump hunting to identify differentially methylated regions
## in epigenetic epidemiology studies. International journal of
## epidemiology, 41(1), 200?209. doi: 10.1093/ije/dyr238.
##
## https://github.com/rafalab/bumphunter
bh.clusterMaker <- function(chr, pos, assumeSorted = FALSE, maxGap=300){
    nonaIndex <- which(!is.na(chr) & !is.na(pos))
    Indexes <- split(nonaIndex, chr[nonaIndex])
    clusterIDs <- rep(NA, length(chr))
    LAST <- 0
    for(i in seq(along = Indexes)){
        Index <- Indexes[[i]]
        x <- pos[Index]
        if(!assumeSorted){
            Index <- Index[order(x)]
            x <- pos[Index]
        }
        y <- as.numeric(diff(x) > maxGap)
        z <- cumsum(c(1, y))
        clusterIDs[Index] <- z + LAST
        LAST <- max(z) + LAST
    }
    clusterIDs
}


## create a variable that is correlated with a CpG site s
## and only correlated with surrounding CpG sites to the
## extent that they are correlated with site s.
generate.fake.bump.var <- function(mat, chr, pos, cluster.sites=40, bump.sites=20, cluster.position=0.5, r=0.5, maxgap=500) {
    stopifnot(bump.sites <= cluster.sites)
    stopifnot(cluster.position >= 0 && cluster.position <= 1)
    stopifnot(r >= 0 && r <= 1)
    ## identify clusters
    clusters <- bh.clusterMaker(chr, pos, maxGap=maxgap)
    cluster.size <- table(clusters)
    
    stopifnot(length(which(cluster.size >= cluster.sites)) > 0)
    
    cluster.size <- cluster.size[which(cluster.size >= cluster.sites)]

    ## select a cluster
    cluster <- sample(names(cluster.size), 1)

    ## indices of cpg sites in the cluster
    cluster.idx <- which(clusters==cluster)

    ## indices of cpg sites in the bump within in the cluster
    bump.center <- max(1, floor(length(cluster.idx)*cluster.position))
    bump.start <- max(1, bump.center - floor(bump.sites/2))
    bump.end <- min(length(cluster.idx), bump.start + bump.sites -1)
    bump.idx <- cluster.idx[bump.start:bump.end]

    ## standardize methylation levels
    mat.scaled <- t(scale(t(mat[bump.idx,])))

    ## identify the middle CpG site, it will be associated with the variable of interest
    site.idx <- which(bump.idx == cluster.idx[bump.center])

    ## calculate correlation of other CpG sites with the middle CpG site
    r.sites <- apply(mat.scaled, 1, function(mat) cor(mat, mat.scaled[site.idx,]))

    ## generate 100 variables associated with the middle CpG site
    vars <- sapply(1:100, function(i) rcor(mat.scaled[site.idx,], r))
    ## calculate correlation of each variable with each cpg site in the 'bump'
    vars.r <- apply(vars, 2, function(var) apply(mat.scaled, 1, function(mat) cor(mat, var)))
    ## pick the variable that most resembles the relationship of the middle CpG site with the other bump cpG sites
    var.idx <- which.max(apply(vars.r, 2, function(var.r) cor(var.r, r.sites)))
    var <- vars[,var.idx]
    list(var=var, bump.idx=bump.idx, cluster.idx=cluster.idx)
}

## create a variable that corresponds to a true bump
generate.true.bump.var <- function(mat, chr, pos, cluster.sites=40, bump.sites=20, cluster.position=0.5, r=0.5, maxgap=500) {
    stopifnot(bump.sites <= cluster.sites)
    stopifnot(cluster.position >= 0 && cluster.position <= 1)
    stopifnot(r >= 0 && r <= 1)
    ## identify clusters
    clusters <- bh.clusterMaker(chr, pos, maxGap=maxgap)
    cluster.size <- table(clusters)
    cluster.size <- cluster.size[which(cluster.size >= cluster.sites)]

    ## select a cluster
    cluster <- sample(names(cluster.size), 1)

    ## indices of cpg sites in the cluster
    cluster.idx <- which(clusters==cluster)

    ## indices of cpg sites in the bump within in the cluster
    bump.center <- max(1, floor(length(cluster.idx)*cluster.position))
    bump.start <- max(1, bump.center - floor(bump.sites/2))
    bump.end <- min(length(cluster.idx), bump.start + bump.sites -1)
    bump.idx <- cluster.idx[bump.start:bump.end]

    ## sum CpG sites after scaling methylation
    bump.mat <- t(scale(t(mat[bump.idx,,drop=F])))
    bump.avg <- colSums(bump.mat)

    ## create a random variable correlated with the sum
    var <- rcor(bump.avg,r)
    
    list(var=var, bump.idx=bump.idx, cluster.idx=cluster.idx)
}

generate.spatial <- function(x, distance) {
    mean.cor <- function(distance) 1/exp(distance/200)
    r <- rnorm(1, mean=mean.cor(distance), sd=0.15)
    if (r < -1) r <- -1
    if (r > 1) r <- 1
    rcor(x, r)
}

generate.methylation <- function(n, pos) {
    cpg.mean <- rep(NA, length(pos))
    cpg.mean[1] <- runif(1)
    cpg.change <- rnorm(length(pos), sd=0.1)
    for (i in 2:length(cpg.mean)) {
        new <- cpg.mean[i-1] + cpg.change[i]
        if (new < 0 | new > 1) cpg.change[i] <- -cpg.change[i]
        cpg.mean[i] <- cpg.mean[i-1] + cpg.change[i]
    }
    
    cpg.sd <- rnorm(length(pos), mean=0.2, sd=0.025)
    
    methylation <- matrix(NA, ncol=n, nrow=length(pos))
    methylation[1,] <- rnorm(ncol(methylation))
    distance <- tail(pos,-1) - head(pos,-1)
    for (i in 2:nrow(methylation))
        methylation[i,] <- generate.spatial(methylation[i-1,], distance[i-1])

    methylation <- t(scale(t(methylation)))
    methylation <- methylation*cpg.sd + cpg.mean
    idx <- which(methylation < 0 | methylation > 1, arr.ind=T)[,1]
    for (i in unique(idx)) {
        med.cpg <- median(methylation[i,])
        min.cpg <- min(methylation[i,])
        max.cpg <- max(methylation[i,])
        factor.min <- med.cpg/(med.cpg - min.cpg)
        factor.max <- (1-med.cpg)/(max.cpg - med.cpg)
        factor <- 1
        if (min.cpg < 0) factor <- min(factor.min, factor)
        if (max.cpg > 1) factor <- min(factor.max, factor)        
        methylation[i,] <- (methylation[i,] - med.cpg)*factor + med.cpg
    }
    methylation
}



