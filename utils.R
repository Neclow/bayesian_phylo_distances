library(phangorn)

read_distance_matrix <- function(filepath) {
  distance_matrix_csv <- read.csv(filepath)
  name <- distance_matrix_csv[, 1]

  distance_matrix <- distance_matrix_csv[2:ncol(distance_matrix_csv)]
  colnames(distance_matrix) <- rownames(distance_matrix) <- name
  distance_matrix <- -as.matrix(distance_matrix)

  key <- order(row.names(distance_matrix))

  distance_matrix_out <- distance_matrix[key, key]

  return(distance_matrix_out)
}

S <- function(times) {
  bf <- c(0.25, 0.25, 0.25, 0.25)
  eig <- edQt(bf = bf)
  P <- eig$vectors %*% diag(exp(eig$values * times)) %*% eig$inv
  P <- P + 1e-12 # Add a small constant to avoid log(0)
  result <- sum(P * log(P) * bf)
  return(result)
}

Sv <- Vectorize(S)

compute_integral <- function(tau, f, F_, num) {
  it <- seq(0.0, tau, length.out = num)
  step_size <- tau / (num - 1)
  h <- numeric(length(it))
  outside <- numeric(length(it))
  h[1] <- F_[1] * S(it[1])
  S_vals <- Sv(it)
  for (t in 2:length(it)) {
    outside[t] <- F_[t] * S(it[t])
    convolution1 <- 0.0
    convolution2 <- 0.0
    j_seq <- 1:(t - 1)
    convolution1 <- sum(h[j_seq] * f[t - j_seq] * step_size)
    convolution2 <- sum(S_vals[j_seq] * f[j_seq] * step_size)
    h[t] <- outside[t] + convolution1 + convolution2
  }
  return(h)
}

entropic <- function(s, theta) {
  return(theta^2 * Sv(s, eig, bf) * exp(-theta * s))
}
