#' Function(s) to perform the clustering over the samples
#'

get_p_clust <- function(family_kl, mu, dis, nc=10, wobs=rep(1,dim(mu)[1]), wsample=rep(1,dim(mu)[2]), cl = NULL) {
    
    # cluster the samples in the latent space if no clustering provided
    if (is.null(cl)) {
        f <- family_kl$linkfun(mu)
        out <- kmeans(t(f), nc, iter.max = 50)
        cl <- out$cluster # cluster indices for each sample
    } else if (typeof(cl)=='list') {
        # old clustering solution provided, so fetch the cluster indices
        if (is.null(cl$cluster))
            stop('argument cl must be a vector of cluster indices or a clustering object returned by k-means.')
        cl <- cl$cluster
    }

    # (re)compute the cluster centers, because they may be different from the ones
    # returned by kmeans if the samples have differing weights
    nc <- max(cl, na.rm=T) # number of clusters (assumes labeling 1,...,nc)
    centers <- matrix(0, nrow=nc, ncol=dim(mu)[1])
    wcluster <- rep(0,nc) # cluster weights
    eps <- 1e-10
    for (j in 1:nc) {
  	    # compute normalized weights within the cluster, 1-eps is for numerical stability
        ind <- which(cl==j)
        ws <- wsample[ind]/sum(wsample[ind])*(1-eps)
    
        # cluster centers and their weights
        centers[j,] <- mu[,ind,drop=F] %*% ws
        wcluster[j] <- sum(wsample[ind]) # unnormalized weight for the jth cluster
    }
    wcluster <- wcluster/sum(wcluster)

    # compute the dispersion parameters for each cluster
    disps <- sapply(1:nc,
                  function(j) {
                  	# compute normalized weights within the cluster, 1-eps is for numerical stability
                    ind <- which(cl == j)
                    ws <- wsample[ind]/sum(wsample[ind])*(1-eps)
                    # dispersion
                    family_kl$discl_fun( mu[,ind,drop=F], dis[ind], wobs, ws )
                  })

    # combine the results
    p <- list(  mu = unname(t(centers)),
                dis = disps,
                weights = wcluster,
                cl = cl)
    return(p)
}
