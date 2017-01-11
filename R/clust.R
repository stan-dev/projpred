#' Function(s) to perform the clustering over the samples
#'

get_p_clust <- function(family_kl, mu, dis, nc=10, wsample=rep(1,dim(mu)[2]), cl = NULL) {

  # TODO
  # THIS FUNCTION WORKS CURRENTLY ONLY FOR GAUSSIAN FAMILY.
  # SHOULD TAKE FAMILY AS AN INPUT AND ACT ACCORDINGLY

  # cluster the mu-samples if no clustering provided
  cl <- cl %ORifNULL% kmeans(t(mu), nc, iter.max = 50)

  if (typeof(cl)=='double') {
    # only cluster-indices provided, so create the list and put them there
    cl <- list(cluster=cl)
  }

  # (re)compute the cluster centers, because they may be different from the ones
  # returned by kmeans if the samples have differing weights
  nc <- max(cl$cluster) # number of clusters (assumes labeling 1,...,nc)
  centers <- matrix(0, nrow=nc, ncol=dim(mu)[1])
  wcluster <- rep(0,nc) # cluster weights
  for (j in 1:nc) {
    ind <- which(cl$cluster==j)
    ws <- wsample[ind]/sum(wsample[ind]) # normalized weights within the cluster
    centers[j,] <- mu[,ind] %*% ws
    wcluster[j] <- sum(wsample[ind]) # unnormalized weight for the jth cluster
  }
  cl$centers <- centers
  wcluster <- wcluster/sum(wcluster)

  # compute the dispersion parameters for each cluster
  disps <- sapply(1:nc,
                  function(cl_ind) {
                    ind <- which(cl$cluster== cl_ind)
                    ws <- wsample[ind]/sum(wsample[ind]) # normalized weights within the cluster
                    if (length(ind) > 1) {
                      mu_mean <- mu[,ind] %*% ws
                      mu_var <- mu[,ind]^2 %*% ws - mu_mean^2
                      sqrt( sum(ws*dis[ind]^2) + mean(mu_var) )
                    } else
                      sqrt(sum(ws*dis[ind]^2))
                  })

  # combine the results
  p <- list(mu = unname(t(cl$centers)),
            dis = disps,
            weights = wcluster)
  list(cl = cl, p = p)
}
