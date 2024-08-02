library(ape)

bme <- function(tree, dists) {
  key <- order(row.names(dists))
  dists <- dists[key, key]
  p_ij <- ape::cophenetic.phylo(tree)
  key <- order(row.names(p_ij))
  p_ij <- p_ij[key, key]
  dp <- dists * 2^(-p_ij)
  return(2 * sum(dp[upper.tri(dp)]))
}

bme_adj <- function(tree, dists, intercept, gradient) {
  key <- order(row.names(dists))
  dists <- dists[key, key]
  p_ij <- ape::cophenetic.phylo(tree)
  key <- order(row.names(p_ij))
  p_ij <- p_ij[key, key]
  dp <- dists * 2^(-p_ij)
  obj <- -2 * sum(dp[upper.tri(dp)])
  return(-1 * (intercept + gradient * obj))
}

run_rax <- function(bl_tr, name) {
  tree_name <- paste0("data/random_trees/tmp_", name, ".tre")
  ape::write.tree(bl_tr, tree_name)
  exec <- paste0(
    "raxml-ng --evaluate --msa data/b10k_100kbp.fasta  --model GTR+G --tree ", tree_name,
    " --brlen scaled --redo  --threads 64 --nofiles"
  )
  output <- invisible(
    system(
      "bash -l",
      input = c("shopt -s expand_aliases", exec),
      intern = TRUE
    )
  )
  matching_index <- grep("Final LogLikelihood:", output, fixed = TRUE)[1]
  log_likelihood_line <- output[matching_index]
  like <- gsub("Final LogLikelihood:", "", log_likelihood_line)
  like <- gsub(" ", "", like)
  like <- as.numeric(like)
  return(like)
}
