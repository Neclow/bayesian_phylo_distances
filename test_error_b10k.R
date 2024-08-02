library(ape)
library(doParallel)
library(GGally)
library(phangorn)
# library(Rsamtools)

source("plot.R")
source("tree.R")
source("utils.R")

aln <- read.FASTA(paste0("data/b10k_100kbp.fasta"))

D <- -as.matrix(dist.ml(phyDat(aln)))
key <- order(row.names(D))
D <- D[key, key]

best_tree <- read.tree("data/main.tre")
raxml_tree <- read.tree("data/b10k_raxml.tre")

fnames <- c("b10k_100", "b10k_60mill", "b10k_600mill_unnorm", "b10k_90mill_unnorm")

# Get all distance matrices
all_dmats <- list()

all_dmats[["b10k_100kbp"]] <- as.vector(D)

all_trees <- list()
all_trees[["best_tree"]] <- best_tree
all_trees[["raxml_tree"]] <- raxml_tree

all_dists_to_best <- list()
all_dists_to_raxml <- list()

for (i in seq_along(fnames)) {
  fname <- fnames[i]

  D_new <- read_distance_matrix(paste("data/", fname, ".csv", sep = ""))

  bme_tree <- fastme.bal(as.dist(-D))

  all_trees[[fname]] <- bme_tree

  all_dists_to_best[[fname]] <- 1 - RF.dist(
    best_tree, bme_tree,
    normalize = TRUE
  )

  all_dists_to_raxml[[fname]] <- 1 - RF.dist(
    raxml_tree, bme_tree,
    normalize = TRUE
  )

  all_dmats[[fname]] <- D_new
}

dmat_df <- as.data.frame(lapply(all_dmats, as.vector))

# Compare distance matrices
ggpairs(
  dmat_df,
  lower = list(continuous = minimal_lmplot),
  diag = list(continuous = minimal_kdeplot)
)

corrcoefs <- cor(dmat_df)

TR <- unroot(all_trees[["b10k_100"]])
TR$edge.length <- rep(1, Nedge(TR))

D_choice <- all_dmats[["b10k_600mill_unnorm"]]

n_trials <- 200
options(digits = 12)

pb <- progress_bar$new(total = n_trials)

out1 <- matrix(NA, nrow = n_trials, ncol = 4)
pb <- progress_bar$new(total = n_trials)
for (i in 1:n_trials) {
  tr <- rNNI(TR, moves = rpois(1, 1) + 1)
  tr <- unroot(tr)
  tr$edge.length <- rep(1, Nedge(tr))
  like1 <- bme(tr, all_dmats[["b10k_100"]])
  like2 <- bme(tr, all_dmats[["b10k_60mill"]])
  like3 <- bme(tr, D_choice)
  like <- run_rax(tr, i)

  pb$tick()

  out1[i, ] <- c(-1 * like1, -1 * like2, -1 * like3, 1 * like)
}

write.table(out1, "trials.Rdata")


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

1 - RF.dist(raxml_tree, mcc_length, normalize = T)



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
  rfd[i] <- 1 - RF.dist(read.tree(text = un_tst[i]), raxml_tree, normalize = TRUE)
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
