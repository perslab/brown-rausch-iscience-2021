compute.sil <- function(x, dist) {
  dist <- as.matrix(dist)

  # Find b
  b.all <- c()
  nearest.clusters.all <- c()

  for (idx.point in 1:length(x)) {
    b <- Inf
    nearest.cluster <- "NA"

    for (cluster in unique(x)) {

      if (cluster != x[idx.point]) {
        idx.cluster <- which(x == cluster)
        b.try <- mean(dist[idx.point, idx.cluster])
        nearest.cluster <- c(nearest.cluster, cluster)[which.min(c(b, b.try))]
        b <- min(b, b.try)

      }
    }
    b.all <- c(b.all, b)
    nearest.clusters.all <- c(nearest.clusters.all, nearest.cluster)
  }

  # Find a
  a.all <- c()
  idx.point <- 1
  for (idx.point in 1:length(x)) {

    cluster <- x[idx.point]
    idx.cluster <- which(x == cluster)[-idx.point]
    a <- mean(dist[idx.point, idx.cluster])
    a.all <- c(a.all, a)
  }

  sil <- (b.all - a.all) / mapply(max, a.all, b.all)
  sil <- cbind(sil, nearest.clusters.all)
  return(sil)

}
