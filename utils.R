read_distance_matrix <- function(filepath) {
  D_new <- read.csv(filepath)
  name <- D_new[, 1]
  D_new <- D_new[2:ncol(D_new)]
  colnames(D_new) <- rownames(D_new) <- name
  D_new <- -as.matrix(D_new)

  key <- order(row.names(D_new))

  D_new <- D_new[key, key]

  return(D_new)
}
