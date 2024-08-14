library(ape)
library(phangorn)

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

bme_neff <- function(tree, dists, neff) {
  key <- order(row.names(dists))
  dists <- dists[key, key]
  p_ij <- cophenetic.phylo(tree)
  key <- order(row.names(p_ij))
  p_ij <- p_ij[key, key]
  dp <- dists * 2^(-p_ij)
  return(-1 * neff * 2 * sum(dp[upper.tri(dp)]))
}

pathnode <- function(phylo, tipsonly = TRUE) {
  di_tr <- dist.nodes(phylo)
  root_tr <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[, 2])][1]

  if (tipsonly == TRUE) {
    roottotippath <- di_tr[
      as.numeric(rownames(di_tr)) == root_tr,
      seq_along(phylo$tip.label)
    ]

    nodesinpath <- sapply(seq_along(phylo$tip.label), function(x) {
      length(Ancestors(phylo, x))
    })
  } else {
    roottotippath <- di_tr[as.numeric(rownames(di_tr)) == root_tr, ]

    nodesinpath <- sapply(
      1:(length(phylo$tip.label) + phylo$Nnode),
      function(x) length(Ancestors(phylo, x))
    )
  }

  return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
}

run_rax <- function(fasta_name, model_name, tree_name, nofiles = FALSE) {
  exec <- paste0(
    "raxml-ng --evaluate",
    " --msa ", fasta_name,
    " --model ", model_name,
    " --tree ", tree_name,
    " --brlen scaled --redo --threads 32"
  )
  if (!nofiles) {
    exec <- paste0(exec, " --nofiles")
  }
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

mcmc <- function(dists, ngen, burnin, thin, intercept, gradient) {
  lp <- numeric(ngen)
  phy <- ape::rmtopology(1, dim(dists)[1], FALSE, tip.label = colnames(dists), br = 1)[[1]]

  current_tree <- phy
  lp[1] <- bme_adj(phy, dists, intercept, gradient)
  naccept <- 0
  k <- 1
  tree_store <- rep("", ngen)
  # pb <- progress_bar$new(
  #   format = "[:bar] :current/:total (:percent) eta: :eta",
  #   total = ngen - 1
  # )

  for (i in 2:ngen) {
    moves <- rpois(1, 1) + 1
    phy_prop <- phangorn::rNNI(current_tree, moves = moves)

    lp_prop <- bme_adj(phy_prop, dists, intercept, gradient)
    ar <- (-lp_prop + lp[i - 1])

    if (log(runif(1)) < ar) {
      naccept <- naccept + 1
      if (i >= burnin && i %% thin == 0) {
        tree_store[k] <- write.tree(phy_prop)
        k <- k + 1
      }
      current_tree <- phy_prop
      lp[i] <- lp_prop
    } else {
      lp[i] <- lp[i - 1]
      if (i >= burnin && i %% thin == 0) {
        tree_store[k] <- write.tree(current_tree)
        k <- k + 1
      }
    }

    # pb$tick()
  }

  return(list(lengths = lp, trees = tree_store[2:k]))
}
