library(ape)
library(dplyr)
library(doParallel)
library(GGally)
library(magrittr)
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
  like3 <- bme(tr, all_dmats[["b10k_600mill_unnorm"]])

  # raxml
  tree_name <- paste0("data/random_trees/tmp_", i, ".tre")
  write.tree(tr, tree_name)
  like <- run_rax(
    fasta_name = "data/b10k_100kbp.fasta",
    model_name = "GTR+G",
    tree_name = tree_name
  )

  pb$tick()

  out1[i, ] <- c(-1 * like1, -1 * like2, -1 * like3, 1 * like)
}

write.table(out1, "results/trials.Rdata")

ggpairs(
  as.data.frame(out1) %>% set_colnames(c("100", "60M", "600M_unnorm", "raxml")),
  lower = list(continuous = minimal_lmplot),
  diag = list(continuous = minimal_kdeplot)
)

lm600mill_coefs <- coef(lm(y ~ x, data = data.frame(x = out1[, 3], y = out1[, 4])))

intercept <- as.numeric(lm600mill_coefs[1])
gradient <- as.numeric(lm600mill_coefs[2])

optimal_tree <- unroot(bme_tree)
optimal_tree$edge.length <- rep(1, Nedge(optimal_tree))
bme_adj(optimal_tree, D_choice, intercept, gradient)

i <- 0
tree_name <- paste0("data/b10k_calibration_trees/", i, ".tre")
write.tree(optimal_tree, tree_name)
like <- run_rax(
  fasta_name = "data/b10k_100kbp.fasta",
  model_name = "GTR+G",
  tree_name = tree_name,
  nofiles = TRUE
)


burnin <- 5000000
registerDoParallel(20)
runs <- 20
ngen <- 20000000
thin <- 2000
count <- 0

for (i in 2:ngen) {
  if (i >= burnin && i %% thin == 0) count <- count + 1
}

out <- foreach(i = 1:runs) %dopar% {
  run <- mcmc(D_choice, ngen, burnin, thin, intercept, gradient)
  return(run)
}

################################################################################

un <- 0
for (i in 1:runs) {
  un <- un + length(unique(out[[i]]$lengths[burnin:ngen]))
}

################################################################################

samp <- seq(1, length(burnin:ngen), length.out = count)
mcmc_matrix <- matrix(ncol = runs, nrow = length(samp))
for (i in 1:runs) {
  mcmc_matrix[, i] <- out[[i]]$lengths[burnin:ngen][samp]
}

trees <- read.tree(text = out[[1]]$trees[1])
for (i in 1:runs) {
  print(i)
  for (j in seq_along(out[[i]]$trees)) {
    trees <- c(trees, read.tree(text = out[[i]]$trees[j]))
  }
}

write.tree(trees, "results/b10k.tree")
write.csv(mcmc_matrix, "results/b10k_mcmc_matrix.csv")

plot(
  mcmc_matrix[, 1],
  type = "l",
  ylim = c(min(mcmc_matrix), max(mcmc_matrix)),
  ylab = "Log Posterior",
  xlab = "Generations"
)

for (i in 2:runs) {
  lines(mcmc_matrix[, i], col = alpha(i, 0.5))
}

mc <- mcc(trees, rooted = FALSE)

mcc_length <- phytools::optim.phylo.ls(
  -as.dist(D_choice),
  stree = mc,
  set.neg.to.zero = TRUE,
  fixed = TRUE,
  tol = 1e-10,
  collapse = FALSE
)

support_tree <- plotBS(mcc_length, trees, type = "phylogram")
write.tree(support_tree, file = "results/b10k_support.tree")

rooted_support_tree <- phytools::midpoint.root(mcc_length)
rooted_support_tree$edge.length[rooted_support_tree$edge.length == 0] <- 1e-4

out_tree <- support_tree
out_tree$node.label <- 100 - out_tree$node.label
write.tree(out_tree, file = "results/out.tree")

print(1 - RF.dist(raxml_tree, mcc_length, normalize = TRUE))

tst <- foreach(i = seq_along(trees), .combine = c) %dopar% {
  tmp <- trees[i]
  tmp <- tmp[[1]]
  tmp$edge.length <- NULL
  tmp <- ladderize(tmp)
  return(write.tree(tmp))
}

un_tst <- unique(tst)
print(length(un_tst))

rfd <- rep(0, length(un_tst))
for (i in seq_along(un_tst)) {
  rfd[i] <- 1 - RF.dist(
    read.tree(text = un_tst[i]), raxml_tree,
    normalize = TRUE
  )
}

nsamp <- 2000
sm <- sample(seq_along(trees), nsamp)
a <- 1:nsamp
b <- 1:nsamp
c <- 1:nsamp
for (i in 1:nsamp) {
  a[i] <- bme(trees[[i]], all_dmats[["b10k_60mill"]])
  b[i] <- bme(trees[[i]], D_choice)
  c[i] <- bme(trees[[i]], all_dmats[["b10k_100"]])
}
