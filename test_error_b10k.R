library(ape)
library(phangorn)
# library(Rsamtools)

bme <- function(tree, dists) {
  key <- order(row.names(dists))
  dists <- dists[key, key]
  p_ij <- ape::cophenetic.phylo(tree)
  key <- order(row.names(p_ij))
  p_ij <- p_ij[key, key]
  dp <- dists * 2^(-p_ij)
  return(2 * sum(dp[upper.tri(dp)]))
}

bme_adj <- function(tree, dists) {
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
    " --brlen scaled --redo  --threads 1 --nofiles"
  )
  output <- invisible(system(exec, intern = TRUE))
  matching_index <- grep("Final LogLikelihood:", output, fixed = TRUE)[1]
  log_likelihood_line <- output[matching_index]
  like <- gsub("Final LogLikelihood:", "", log_likelihood_line)
  like <- gsub(" ", "", like)
  like <- as.numeric(like)
  return(like)
}

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

aln <- read.FASTA(paste0("data/b10k_100kbp.fasta"))

chars <- as.character(aln)

D <- -as.matrix(dist.ml(phyDat(aln)))
key <- order(row.names(D))
D <- D[key, key]

best_tree <- read.tree("data/main.tre")
raxml <- read.tree("data/b10k_raxml.tre")

D_new <- read_distance_matrix("data/b10k_100.csv")

D_new2 <- read_distance_matrix("data/b10k_60mill.csv")

D_new3 <- read_distance_matrix("data/b10k_600mill_unnorm.csv")

D_new4 <- read_distance_matrix("data/b10k_90mill_unnorm.csv")

bme_tree <- fastme.bal(as.dist(-D))
bme_tree2 <- fastme.bal(as.dist(-D_new))
bme_tree3 <- fastme.bal(as.dist(-D_new2))
bme_tree4 <- fastme.bal(as.dist(-D_new3))
bme_tree5 <- fastme.bal(as.dist(-D_new4))

1 - RF.dist(best_tree, bme_tree, normalize = TRUE)
1 - RF.dist(best_tree, bme_tree2, normalize = TRUE)
1 - RF.dist(best_tree, bme_tree3, normalize = TRUE)
1 - RF.dist(best_tree, bme_tree4, normalize = TRUE)
1 - RF.dist(best_tree, bme_tree5, normalize = TRUE)

1 - RF.dist(raxml, best_tree, normalize = TRUE)
1 - RF.dist(raxml, bme_tree, normalize = TRUE)
1 - RF.dist(raxml, bme_tree2, normalize = TRUE)
1 - RF.dist(raxml, bme_tree3, normalize = TRUE)
1 - RF.dist(raxml, bme_tree4, normalize = TRUE)
1 - RF.dist(raxml, bme_tree5, normalize = TRUE)

par(mfrow = c(1, 4))
plot(as.vector(D_new), as.vector(D_new2))
abline(0, 1, col = "red")
lm(x ~ y, data.frame(x = as.vector(D_new), y = as.vector(D_new2)))


plot(as.vector(D_new), as.vector(D_new3))
abline(0, 1, col = "red")
lm(x ~ y, data.frame(x = as.vector(D_new), y = as.vector(D_new3)))
lm_adjust <- lm(x ~ y, data.frame(x = as.vector(D_new), y = as.vector(D_new3)))


plot(as.vector(D_new2), as.vector(D_new3))
abline(0, 1, col = "red")
lm(x ~ y, data.frame(x = as.vector(D_new2), y = as.vector(D_new3)))


plot(as.vector(D_new), as.vector(D_new4))
abline(0, 1, col = "red")
lm(x ~ y, data.frame(x = as.vector(D_new), y = as.vector(D_new4)))


TR <- bme_tree2
TR <- unroot(TR)
TR$edge.length <- rep(1, Nedge(TR))


D_choice <- D_new3

ntrials <- 200
options(digits = 12)
setwd("~/tools/raxml-ng_v1.2.0_linux_x86_64/")
out1 <- matrix(NA, nrow = ntrials, ncol = 4)
library(doParallel)
registerDoParallel(50)
out1 <- foreach(i = 1:ntrials, .combine = rbind) %dopar% {
  tr <- rNNI(TR, moves = rpois(1, 1) + 1)
  tr <- unroot(tr)
  tr$edge.length <- rep(1, Nedge(tr))
  like1 <- bme(tr, D_new)
  like2 <- bme(tr, D_new2)
  like3 <- bme(tr, D_choice)
  like <- run_rax(tr, i)
  return(cbind(-1 * like1, -1 * like2, -1 * like3, like))
}


par(mfrow = c(1, 3))

lm100k <- lm(y ~ x, data = data.frame(x = out1[, 1], y = out1[, 4]))
plot(coef(lm100k)[1] + coef(lm100k)[2] * out1[, 1], out1[, 4], pch = 16)
abline(0, 1, col = "red")

lm60mill <- lm(y ~ x, data = data.frame(x = out1[, 2], y = out1[, 4]))
plot(coef(lm60mill)[1] + coef(lm60mill)[2] * out1[, 2], out1[, 4], pch = 16)
abline(0, 1, col = "red")

lm600mill <- lm(y ~ x, data = data.frame(x = out1[, 3], y = out1[, 4]))
plot(coef(lm600mill)[1] + coef(lm600mill)[2] * out1[, 3], out1[, 4], pch = 16)
abline(0, 1, col = "red")

lm100k
lm60mill
lm600mill

intercept <- as.numeric(coef(lm600mill)[1])
gradient <- as.numeric(coef(lm600mill)[2])

optimal_tree <- bme_tree
optimal_tree <- unroot(optimal_tree)
optimal_tree$edge.length <- rep(1, Nedge(optimal_tree))
bme_adj(optimal_tree, D_choice)


i <- 0
write.tree(optimal_tree, paste0("data/b10k_calibration_trees/", i, ".tre"))
exec <- paste0("raxml-ng --evaluate --msa data/b10k_100kbp.fasta --model GTR+G --tree data/b10k_calibration_trees/", i, ".tre --brlen scaled --redo")
invisible(system(exec, intern = TRUE))
# 	system(exec)
log_file_path <- "data/b10k_100kbp.fasta.raxml.log"
grep_command <- paste("grep 'Final LogLikelihood:'", log_file_path, sep = " ")
log_likelihood_line <- system(grep_command, intern = TRUE)
like <- gsub("Final LogLikelihood:", "", log_likelihood_line)
like <- gsub(" ", "", like)
like <- as.numeric(like)
like


burnin <- 5000000
library(doParallel)
registerDoParallel(20)
runs <- 20
ngen <- 20000000
thin <- 2000
count <- 0
for (i in 2:ngen) {
  if (i >= burnin & i %% thin == 0) count <- count + 1
}
count * 20


mcmc <- function(dists, ngen) {
  lp <- numeric(ngen)
  trees <- vector("list", ngen)
  phy <- ape::rmtopology(1, dim(dists)[1], FALSE, tip.label = colnames(dists), br = 1)[[1]]

  current_tree <- phy
  lp[1] <- bme_adj(phy, dists)
  naccept <- 0
  k <- 1
  tree_store <- rep("", ngen)
  for (i in 2:ngen) {
    moves <- rpois(1, 1) + 1
    phy_prop <- phangorn::rNNI(current_tree, moves = moves)

    lp_prop <- bme_adj(phy_prop, dists)
    AR <- (-lp_prop + lp[i - 1])

    if (log(runif(1)) < AR) {
      naccept <- naccept + 1
      if (i >= burnin & i %% thin == 0) {
        tree_store[k] <- write.tree(phy_prop)
        k <- k + 1
      }
      current_tree <- phy_prop
      lp[i] <- lp_prop
    } else {
      lp[i] <- lp[i - 1]
      if (i >= burnin & i %% thin == 0) {
        tree_store[k] <- write.tree(current_tree)
        k <- k + 1
      }
    }
  }

  return(list(lengths = lp, trees = tree_store[2:k]))
}

out <- foreach(i = 1:runs) %dopar% {
  run <- mcmc(D_choice, ngen)
  return(run)
}


################################################################################

un <- 0
for (i in 1:runs) {
  un <- un + length(unique(out[[i]]$lengths[burnin:ngen]))
}

################################################################################

samp <- seq(1, length(burnin:ngen), length.out = count)
mcmc_matrix <- matrix(nrow = runs, ncol = length(samp))
for (i in 1:runs) {
  mcmc_matrix[i, ] <- out[[i]]$lengths[burnin:ngen][samp]
}
length(unique(as.vector(mcmc_matrix)))

library(scales)
plot(mcmc_matrix[1, ], type = "l", ylim = c(min(mcmc_matrix), max(mcmc_matrix)), ylab = "Log Posterior", xlab = "Generations")
for (i in 2:runs) {
  lines(mcmc_matrix[i, ], col = alpha(i, 0.5))
}


trees <- read.tree(text = out[[1]]$trees[1])
for (i in 1:runs) {
  print(i)
  for (j in 1:length(out[[i]]$trees)) {
    trees <- c(trees, read.tree(text = out[[i]]$trees[j]))
  }
}
mc <- mcc(trees, rooted = F)

mcc_length <- phytools::optim.phylo.ls(-1 * D_choice, stree = mc, set.neg.to.zero = TRUE, fixed = TRUE, tol = 1e-10, collapse = FALSE)

support_tree <- plotBS(mcc_length, trees, type = "phylogram")
write.tree(support_tree)

rooted_support_tree <- phytools::midpoint.root(mcc_length)
rooted_support_tree$edge.length[rooted_support_tree$edge.length == 0] <- 1e-4

out_tree <- support_tree
out_tree$node.label <- 100 - out_tree$node.label
write.tree(out_tree)

1 - RF.dist(raxml, mcc_length, normalize = T)



write.tree(trees, paste0("data/b10k.tree"))
write.csv(mcmc_matrix, paste0("data/b10k_mcmc_matrix.csv"))

tst <- foreach(i = 1:length(trees), .combine = c) %dopar% {
  tmp <- trees[i]
  tmp <- tmp[[1]]
  tmp$edge.length <- NULL
  tmp <- ladderize(tmp)
  return(write.tree(tmp))
}

un_tst <- unique(tst)
length(un_tst)


rfd <- rep(0, length(un_tst))
for (i in 1:length(un_tst)) {
  rfd[i] <- 1 - RF.dist(read.tree(text = un_tst[i]), raxml, normalize = TRUE)
}


nsamp <- 2000
sm <- sample(1:length(trees), nsamp)
a <- 1:nsamp
b <- 1:nsamp
c <- 1:nsamp
for (i in 1:nsamp) {
  a[i] <- bme(trees[[i]], D_new2)
  b[i] <- bme(trees[[i]], D_choice)
  c[i] <- bme(trees[[i]], D_new)
}
