library(ape)
# library(bayesm)
library(doParallel)
library(ggplot2)
library(gridExtra)
# library(fitdistrplus)
library(NELSI)
library(phangorn)
library(TreeSim)

source("tree.R")
source("utils.R")

ntax <- 50
neff <- 5000

registerDoParallel(60)
count <- 1
slow <- TRUE

rates <- c(
  seq(1, 2.5, length = 4),
  seq(0.1, 0.9, length = 9),
  seq(0.025, 0.075, length = 3),
  seq(0.005, 0.01, length = 3)
)

probs <- c(0.025, 0.5, 0.975)

a1 <- a2 <- a3 <- a4 <- a5 <- a6 <- a7 <- a8 <- a9 <- matrix(
  NA,
  nrow = length(rates),
  ncol = 3
)

get_tree_and_dmat_from_sim <- function(rate) {
  meanrtt <- mean(pathnode(realtree)[[1]])
  realtree$edge.length <- (realtree$edge.length * rate) / meanrtt
  aln <- simSeq(realtree, l = neff)
  write.FASTA(as.DNAbin(aln), paste0("tmp.fa"))

  D_orig <- dist.ml(aln, model = "JC69", exclude = "pairwise")

  bme_tree <- fastme.bal(D_orig)

  D_orig <- as.matrix(D_orig)
  key <- order(row.names(D_orig))
  D_orig <- D_orig[key, key]

  return(list(aln = aln, D_orig = D_orig, bme_tree = bme_tree))
}

make_D_new <- function(D_orig, dsvals, it) {
  D_new <- matrix(0, nrow = nrow(D_orig), ncol = ncol(D_orig))
  for (n in seq_len(nrow(D_orig))) {
    for (m in seq_len(ncol(D_orig))) {
      if (n != m) {
        wh <- abs(it - D_orig[n, m])
        D_new[n, m] <- dsvals[which.min(wh)]
      }
    }
  }
  rownames(D_new) <- rownames(D_orig)
  colnames(D_new) <- colnames(D_orig)

  return(D_new)
}

fit_trees <- function(realtree_unrooted, aln, D_orig, D_new, i) {
  tr <- rSPR(realtree_unrooted, moves = rpois(1, 1) + 1)

  tr$edge.length <- rep(1, Nedge(tr))
  tr <- unroot(tr)
  theta <- 2 * (ntax - 2) / bme(tr, D_orig)

  fit_jc <- pml(tr, aln, control = pml.control(trace = 0))
  fit_jc <- optim.pml(
    fit_jc,
    optNni = FALSE,
    optBf = FALSE,
    optEdge = TRUE,
    control = pml.control(trace = 0)
  )

  if (slow) {
    D_step <- D_orig

    maxim <- round(max(D_step), 2) + 0.05
    num <- 2000
    it <- seq(0.0, maxim, length.out = num)
    f <- dexp(it, theta)
    F_ <- 1 - pexp(it, theta)

    dsvals <- compute_integral(maxim, f, F_, num)

    D_mid <- make_D_new(D_orig = D_step, dsvals = dsvals, it = it)
  } else {
    b <- integrate(entropic, 0, Inf, theta)$value
    D_mid <- b * D_orig
  }

  like1 <- bme_neff(tr, D_mid, neff)
  like2 <- bme_neff(tr, D_new, neff)

  tree_name <- paste0("data/random_trees/tmp_", i, ".tre")
  write.tree(tr, tree_name)

  like1r <- run_rax(
    fasta_name = "tmp.fa",
    model_name = "JC",
    tree_name = tree_name
  )

  fit_jc <- pml(tr, aln, control = pml.control(trace = 0))
  fit_jc <- optim.pml(
    fit_jc,
    optNni = FALSE,
    optBf = FALSE,
    optEdge = TRUE,
    control = pml.control(trace = 0)
  )

  return(cbind(-1 * like1, like1r, -1 * like2))
}

total_replicates <- 100

for (rate in rates) {
  tmp <- matrix(NA, nrow = total_replicates, ncol = 9)
  for (zzz in 1:total_replicates) {
    chronogram <- ladderize(sim.bd.taxa.age(ntax, 1, 0.5, 0.1, 1, 65)[[1]])
    meanrate <- 0.05
    ratevar <- 0.05
    phylogram <- simulate.white.noise(
      chronogram,
      params = list(rate = meanrate, var = ratevar)
    )[[1]]
    realtree <- phylogram

    sim_data <- get_tree_and_dmat_from_sim(rate)

    aln <- sim_data[["aln"]]
    D_orig <- sim_data[["D_orig"]]
    bme_tree <- sim_data[["bme_tree"]]

    tr <- bme_tree
    tr$edge.length <- rep(1, Nedge(tr))
    theta <- 2 * (ntax - 2) / bme(tr, D_orig)

    maxim <- round(max(D_orig), 2) + 0.01
    num <- 2000
    it <- seq(0.0, maxim, length.out = num)
    f <- dexp(it, theta)
    F_ <- 1 - pexp(it, theta)

    dsvals <- compute_integral(maxim, f, F_, num)

    D_new <- make_D_new(D_orig, dsvals, it)

    bme_tree2 <- fastme.bal(as.dist(-1 * D_new))
    dist_bme <- 1 - RF.dist(bme_tree, realtree, normalize = TRUE)
    dist_entropic <- 1 - RF.dist(bme_tree2, realtree, normalize = TRUE)

    fit_jc <- pml(bme_tree, aln, control = pml.control(trace = 0))
    fit_jc <- optim.pml(
      fit_jc,
      optNni = TRUE,
      optBf = FALSE,
      optEdge = TRUE,
      control = pml.control(trace = 0)
    )

    dist_like <- 1 - RF.dist(fit_jc$tree, realtree, normalize = TRUE)

    bme_tree2$edge.length <- rep(1, Nedge(bme_tree2))

    realtree_unrooted <- realtree
    realtree_unrooted$edge.length <- rep(1, Nedge(realtree_unrooted))
    realtree_unrooted <- unroot(realtree_unrooted)
    bme_neff(realtree_unrooted, D_new, neff)

    ntrials <- total_replicates
    out1 <- foreach(i = 1:ntrials, .combine = rbind) %dopar% {
      fit_trees(realtree_unrooted, aln, D_orig, D_new, i)
    }

    lm_model <- lm(y ~ x, data = data.frame(x = out1[, 1], y = out1[, 2]))
    lm_model2 <- lm(y ~ x, data = data.frame(x = out1[, 3], y = out1[, 2]))

    mae1 <- mean(
      abs(
        coef(lm_model)[1] + coef(lm_model)[2] * out1[, 1] - out1[, 2]
      ) / -out1[, 2]
    )
    mae2 <- mean(
      abs(
        coef(lm_model2)[1] + coef(lm_model2)[2] * out1[, 3] - out1[, 2]
      ) / -out1[, 2]
    )

    sd1 <- sd(
      abs(
        coef(lm_model)[1] + coef(lm_model)[2] * out1[, 1] - out1[, 2]
      ) / -out1[, 2]
    )
    sd2 <- sd(
      abs(
        coef(lm_model2)[1] + coef(lm_model2)[2] * out1[, 3] - out1[, 2]
      ) /
        -out1[, 2]
    )

    tmp[zzz, ] <- c(
      coef(lm_model)[2],
      coef(lm_model2)[2],
      dist_like, dist_bme, dist_entropic,
      mae1, mae2,
      sd1, sd2
    )
  }
  a1[count, ] <- quantile(tmp[, 1], probs = probs)
  a2[count, ] <- quantile(tmp[, 2], probs = probs)
  a3[count, ] <- quantile(tmp[, 3], probs = probs)
  a4[count, ] <- quantile(tmp[, 4], probs = probs)
  a5[count, ] <- quantile(tmp[, 5], probs = probs)
  a6[count, ] <- quantile(tmp[, 6], probs = probs)
  a7[count, ] <- quantile(tmp[, 7], probs = probs)
  a8[count, ] <- quantile(tmp[, 8], probs = probs)
  a9[count, ] <- quantile(tmp[, 9], probs = probs)
  count <- count + 1
  print(count)
}

# Your data
coeff_entropic <- data.frame(
  x = rates,
  median = a1[, 2],
  low_ci = a1[, 1],
  high_ci = a1[, 3]
)
coeff_bme <- data.frame(
  x = rates,
  median = a2[, 2],
  low_ci = a2[, 1],
  high_ci = a2[, 3]
)

rf_like <- data.frame(
  x = rates,
  median = a3[, 2],
  low_ci = a3[, 1],
  high_ci = a3[, 3]
)

rf_bme <- data.frame(
  x = rates,
  median = a4[, 2],
  low_ci = a4[, 1],
  high_ci = a4[, 3]
)

rf_entropic <- data.frame(
  x = rates,
  median = a5[, 2],
  low_ci = a5[, 1],
  high_ci = a5[, 3]
)

mae_entropic <- data.frame(
  x = rates,
  median = a6[, 2],
  low_ci = a6[, 1],
  high_ci = a6[, 3]
)

mae_bme <- data.frame(
  x = rates,
  median = a7[, 2],
  low_ci = a7[, 1],
  high_ci = a7[, 3]
)

sd_entropic <- data.frame(
  x = rates,
  median = a8[, 2],
  low_ci = a8[, 1],
  high_ci = a8[, 3]
)

sd_bme <- data.frame(
  x = rates,
  median = a9[, 2],
  low_ci = a9[, 1],
  high_ci = a9[, 3]
)

# Plotting with ggplot2
p1 <- ggplot() +
  geom_point(data = coeff_entropic, aes(x = x, y = median), color = "blue", alpha = 0.5, size = 3) +
  geom_errorbar(
    data = coeff_entropic,
    aes(
      ymin = low_ci,
      ymax = high_ci,
      y = median,
      x = x
    ),
    width = 0.05,
    color = "blue",
    linewidth = 1,
    alpha = 0.25
  ) +
  geom_point(
    data = coeff_bme,
    aes(x = x, y = median),
    color = "red",
    alpha = 0.5,
    size = 3
  ) +
  geom_errorbar(
    data = coeff_bme,
    aes(
      ymin = low_ci,
      ymax = high_ci,
      y = median,
      x = x
    ),
    width = 0.05,
    color = "red",
    linewidth = 1,
    alpha = 0.25
  ) +
  geom_vline(aes(xintercept = 0.1), linetype = "dotted", alpha = 0.75) +
  geom_vline(aes(xintercept = 0.5), linetype = "dotted", alpha = 0.75) +
  labs(x = "Rate", title = "Gradient", y = "g") +
  theme_bw() +
  scale_x_log10() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12), # Adjust the size for tick labels
    plot.title = element_text(size = 16)
  )

p2 <- ggplot() +
  geom_point(data = rf_like, aes(x = x, y = median), color = "black", alpha = 0.5, size = 3) +
  geom_errorbar(
    data = rf_like,
    aes(ymin = low_ci, ymax = high_ci, y = median, x = x),
    width = 0.075,
    color = "black",
    size = 1,
    alpha = 0.25
  ) +
  geom_point(
    data = rf_entropic,
    aes(x = x, y = median),
    color = "blue",
    alpha = 0.5,
    size = 3
  ) +
  geom_errorbar(
    data = rf_entropic,
    aes(ymin = low_ci, ymax = high_ci, y = median, x = x),
    width = 0.05,
    color = "blue",
    size = 1,
    alpha = 0.25
  ) +
  geom_point(
    data = rf_bme,
    aes(x = x, y = median),
    color = "red",
    alpha = 0.5,
    size = 3
  ) +
  geom_errorbar(
    data = rf_bme,
    aes(ymin = low_ci, ymax = high_ci, y = median, x = x),
    width = 0.05,
    color = "red",
    size = 1,
    alpha = 0.25
  ) +
  geom_vline(aes(xintercept = 0.1), linetype = "dotted", alpha = 0.75) +
  geom_vline(aes(xintercept = 0.5), linetype = "dotted", alpha = 0.75) +
  labs(
    x = "Rate",
    y = "1 - Normalised RF distance",
    title = "Topological accuracy"
  ) +
  theme_bw() +
  scale_x_log10() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12), # Adjust the size for tick labels
    plot.title = element_text(size = 16)
  )

p3 <- ggplot() +
  geom_point(
    data = mae_entropic,
    aes(x = x, y = 100 * median),
    color = "blue",
    alpha = 0.5,
    size = 3
  ) +
  geom_errorbar(
    data = mae_entropic,
    aes(ymin = 100 * low_ci, ymax = 100 * high_ci, y = median, x = x),
    width = 0.05,
    color = "blue",
    size = 1,
    alpha = 0.25
  ) +
  geom_point(
    data = mae_bme,
    aes(x = x, y = 100 * median),
    color = "red",
    alpha = 0.5,
    size = 3
  ) +
  geom_errorbar(
    data = mae_bme,
    aes(ymin = 100 * low_ci, ymax = 100 * high_ci, y = median, x = x),
    width = 0.05,
    color = "red",
    size = 1,
    alpha = 0.25
  ) +
  geom_vline(aes(xintercept = 0.1), linetype = "dotted", alpha = 0.75) +
  geom_vline(aes(xintercept = 0.5), linetype = "dotted", alpha = 0.75) +
  labs(x = "Rate", y = "MAPE (%)", title = "Entropic vs Felsenstein log like") +
  theme_bw() +
  scale_x_log10() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12), # Adjust the size for tick labels
    plot.title = element_text(size = 16)
  )

grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
