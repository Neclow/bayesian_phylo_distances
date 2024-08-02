library(phangorn)
library(fitdistrplus)

library(phangorn)
library(bayesm)
library(TreeSim)
library(fitdistrplus)
library(NELSI)

compute_integral <- function(tau, f, F, num) {
  it <- seq(0.0, tau, length.out = num)
  step_size <- tau / (num - 1)
  h <- numeric(length(it))
  int1 <- numeric(length(it))
  int2 <- numeric(length(it))
  outside <- numeric(length(it))
  h[1] <- F[1] * S(it[1])
  S_vals <- Sv(it)
  for (t in 2:length(it)) {
    outside[t] <- F[t] * S(it[t])
    convolution1 <- 0.0
    convolution2 <- 0.0
    j_seq <- 1:(t - 1)
    convolution1 <- sum(h[j_seq] * f[t - j_seq] * step_size)
    convolution2 <- sum(S_vals[j_seq] * f[j_seq] * step_size)
    h[t] <- outside[t] + convolution1 + convolution2
  }
  return(h)
}
S <- function(times) {
  P <- eig$vectors %*% diag(exp(eig$values * times)) %*% eig$inv
  P <- P + 1e-12 # Add a small constant to avoid log(0)
  log_P <- log(P)
  result <- sum(P * log_P * bf)
  return(result)
}
Sv <- Vectorize(S)
entropic <- function(s, theta) {
  return(theta^2 * Sv(s) * exp(-theta * s))
}

bme <- function(tree, dists) {
  key <- order(row.names(dists))
  dists <- dists[key, key]
  p_ij <- ape::cophenetic.phylo(tree)
  key <- order(row.names(p_ij))
  p_ij <- p_ij[key, key]
  dp <- dists * 2^(-p_ij)
  return(-1 * neff * 2 * sum(dp[upper.tri(dp)]))
}
bme_orig <- function(tree, dists) {
  key <- order(row.names(dists))
  dists <- dists[key, key]
  p_ij <- ape::cophenetic.phylo(tree)
  key <- order(row.names(p_ij))
  p_ij <- p_ij[key, key]
  dp <- dists * 2^(-p_ij)
  return(2 * sum(dp[upper.tri(dp)]))
}

pathnode <- function(phylo, tipsonly = T) {
  require(phangorn)
  di.tr <- dist.nodes(phylo)
  root.tr <- phylo$edge[, 1][!(phylo$edge[, 1] %in% phylo$edge[
    ,
    2
  ])][1]
  tr.depth <- max(di.tr[as.numeric(colnames(di.tr)) == root.tr, ])
  if (tipsonly == TRUE) {
    roottotippath <- di.tr[as.numeric(rownames(di.tr)) ==
      root.tr, 1:length(phylo$tip.label)]
    nodesinpath <- sapply(1:length(phylo$tip.label), function(x) {
      length(Ancestors(
        phylo,
        x
      ))
    })
  } else {
    roottotippath <- di.tr[as.numeric(rownames(di.tr)) ==
      root.tr, ]
    nodesinpath <- sapply(
      1:(length(phylo$tip.label) + phylo$Nnode),
      function(x) length(Ancestors(phylo, x))
    )
  }
  # plot(roottotippath, nodesinpath, xlab = "Root-to-tip path length", ylab = "Number of parent nodes", pch = 20)
  return(list(roottotippath = roottotippath, nodesinpath = nodesinpath))
}

bf <- c(0.25, 0.25, 0.25, 0.25)
eig <- edQt(bf = bf)

ntax <- 50
neff <- 5000



library(doParallel)
registerDoParallel(60)
count <- 1
SLOW <- TRUE
rates <- c(2.5, 2, 1.5, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.005)
a1 <- a2 <- a3 <- a4 <- a5 <- a6 <- a7 <- a8 <- a9 <- matrix(NA, nrow = length(rates), ncol = 3)

total_replicates <- 100

run_iqtree <- function(bl_tr) {
  write.tree(bl_tr, paste0("/home/sjbhatt/tmp.tre"))
  exec <- "iqtree2 -s /home/sjbhatt/tmp.fa -z /home/sjbhatt/tmp.tre -m JC -n 0 -redo"
  invisible(system(exec, intern = TRUE))
  log_file_path <- "/home/sjbhatt/tmp.fa.log"
  grep_command <- paste("grep 'Optimal log-likelihood:'", log_file_path, sep = " ")
  log_likelihood_line <- system(grep_command, intern = TRUE)
  like <- gsub("Optimal log-likelihood:", "", log_likelihood_line)
  like <- gsub(" ", "", like)
  like <- as.numeric(like)
  return(like)
}


run_rax <- function(bl_tr, name) {
  tree_name <- paste0("/home/sjbhatt/random_trees/tmp_", name, ".tre")
  write.tree(bl_tr, tree_name)
  exec <- paste0("/home/sjbhatt/raxml-ng_v1.2.0_linux_x86_64/./raxml-ng --evaluate --msa /home/sjbhatt/tmp.fa --model JC --tree ", tree_name, " --brlen scaled --opt-model off --opt-branches on --redo  --threads 1 --nofiles")
  output <- invisible(system(exec, intern = TRUE))
  matching_index <- grep("Final LogLikelihood:", output, fixed = TRUE)[1]
  log_likelihood_line <- output[matching_index]
  like <- gsub("Final LogLikelihood:", "", log_likelihood_line)
  like <- gsub(" ", "", like)
  like <- as.numeric(like)
  return(like)
}



for (rate in rates) {
  tmp <- matrix(NA, nrow = total_replicates, ncol = 9)
  for (zzz in 1:total_replicates) {
    # 	realtree=rtree(ntax,br=rexp)

    chronogram <- ladderize(sim.bd.taxa.age(ntax, 1, 0.5, 0.1, 1, 65)[[1]])
    meanrate <- 0.05
    ratevar <- 0.05
    phylogram <- simulate.white.noise(chronogram, params = list(rate = meanrate, var = ratevar))[[1]]
    realtree <- phylogram

    meanrtt <- mean(pathnode(realtree)[[1]])
    realtree$edge.length <- (realtree$edge.length * rate) / meanrtt
    aln <- simSeq(realtree, l = neff)
    write.FASTA(as.DNAbin(aln), paste0("/home/sjbhatt/tmp.fa"))
    chars <- as.character(aln)
    bme_tree <- fastme.bal(dist.ml(aln, model = "JC69", exclude = "pairwise"))
    D_orig <- as.matrix(dist.ml(aln, model = "JC69", exclude = "pairwise"))
    key <- order(row.names(D_orig))
    D_orig <- D_orig[key, key]

    tr <- bme_tree
    tr$edge.length <- rep(1, Nedge(tr))
    theta <- 2 * (ntax - 2) / (bme_orig(tr, D_orig))

    maxim <- round(max(D_orig), 2) + 0.01
    num <- 2000
    it <- seq(0.0, maxim, length.out = num)
    f <- dexp(it, theta)
    F <- 1 - pexp(it, theta)

    dsvals <- compute_integral(maxim, f, F, num)


    D_new <- matrix(0, nrow = nrow(D_orig), ncol = ncol(D_orig))
    for (n in 1:nrow(D_orig)) {
      for (m in 1:ncol(D_orig)) {
        if (n != m) {
          wh <- abs(it - D_orig[n, m])
          D_new[n, m] <- dsvals[which.min(wh)]
        }
      }
    }
    rownames(D_new) <- rownames(D_orig)
    colnames(D_new) <- colnames(D_orig)


    bme_tree2 <- fastme.bal(as.dist(-1 * D_new))
    dist_bme <- 1 - RF.dist(bme_tree, realtree, normalize = T)
    dist_entropic <- 1 - RF.dist(bme_tree2, realtree, normalize = T)

    fitJC <- pml(bme_tree, aln, control = pml.control(trace = 0))
    fitJC <- optim.pml(fitJC, optNni = TRUE, optBf = FALSE, optEdge = TRUE, control = pml.control(trace = 0))

    dist_like <- 1 - RF.dist(fitJC$tree, realtree, normalize = T)

    bme_tree2$edge.length <- rep(1, Nedge(bme_tree2))

    bme <- function(tree, dists) {
      key <- order(row.names(dists))
      dists <- dists[key, key]
      p_ij <- ape::cophenetic.phylo(tree)
      key <- order(row.names(p_ij))
      p_ij <- p_ij[key, key]
      dp <- dists * 2^(-p_ij)
      return(-1 * neff * 2 * sum(dp[upper.tri(dp)]))
    }


    TR <- realtree
    TR$edge.length <- rep(1, Nedge(TR))
    TR <- unroot(TR)
    bme(TR, D_new)

    ntrials <- total_replicates
    out1 <- foreach(i = 1:ntrials, .combine = rbind) %dopar% {
      tr <- rSPR(TR, moves = rpois(1, 1) + 1)

      tr$edge.length <- rep(1, Nedge(tr))
      tr <- unroot(tr)
      theta <- 2 * (ntax - 2) / (bme_orig(tr, D_orig))

      fitJC <- pml(tr, aln, control = pml.control(trace = 0))
      fitJC <- optim.pml(fitJC, optNni = FALSE, optBf = FALSE, optEdge = TRUE, control = pml.control(trace = 0))

      if (SLOW) {
        D_step <- D_orig

        maxim <- round(max(D_step), 2) + 0.05
        num <- 2000
        it <- seq(0.0, maxim, length.out = num)
        f <- dexp(it, theta)
        F <- 1 - pexp(it, theta)

        dsvals <- compute_integral(maxim, f, F, num)

        D_mid <- matrix(0, nrow = nrow(D_step), ncol = ncol(D_step))
        for (n in 1:nrow(D_step)) {
          for (m in 1:ncol(D_step)) {
            if (n != m) {
              wh <- abs(it - D_step[n, m])
              D_mid[n, m] <- dsvals[which.min(wh)]
            }
          }
        }
        rownames(D_mid) <- rownames(D_step)
        colnames(D_mid) <- colnames(D_step)
      } else {
        b <- integrate(entropic, 0, Inf, theta)$value
        D_mid <- b * D_orig
      }

      like1 <- bme(tr, D_mid)
      like2 <- bme(tr, D_new)


      like1r <- run_rax(tr, i)
      fitJC <- pml(tr, aln, control = pml.control(trace = 0))
      fitJC <- optim.pml(fitJC, optNni = FALSE, optBf = FALSE, optEdge = TRUE, control = pml.control(trace = 0))
      like2r <- fitJC$logLik


      return(cbind(-1 * like1, like1r, -1 * like2))
    }

    lm.model <- lm(y ~ x, data = data.frame(x = out1[, 1], y = out1[, 2]))
    lm.model2 <- lm(y ~ x, data = data.frame(x = out1[, 3], y = out1[, 2]))

    mae1 <- mean(abs((coef(lm.model)[1] + coef(lm.model)[2] * out1[, 1]) - out1[, 2]) / -out1[, 2])
    mae2 <- mean(abs((coef(lm.model2)[1] + coef(lm.model2)[2] * out1[, 3]) - out1[, 2]) / -out1[, 2])

    sd1 <- sd(abs((coef(lm.model)[1] + coef(lm.model)[2] * out1[, 1]) - out1[, 2]) / -out1[, 2])
    sd2 <- sd(abs((coef(lm.model2)[1] + coef(lm.model2)[2] * out1[, 3]) - out1[, 2]) / -out1[, 2])


    tmp[zzz, ] <- c(coef(lm.model)[2], coef(lm.model2)[2], dist_like, dist_bme, dist_entropic, mae1, mae2, sd1, sd2)
  }
  a1[count, ] <- quantile(tmp[, 1], probs = c(0.025, 0.5, 0.975))
  a2[count, ] <- quantile(tmp[, 2], probs = c(0.025, 0.5, 0.975))
  a3[count, ] <- quantile(tmp[, 3], probs = c(0.025, 0.5, 0.975))
  a4[count, ] <- quantile(tmp[, 4], probs = c(0.025, 0.5, 0.975))
  a5[count, ] <- quantile(tmp[, 5], probs = c(0.025, 0.5, 0.975))
  a6[count, ] <- quantile(tmp[, 6], probs = c(0.025, 0.5, 0.975))
  a7[count, ] <- quantile(tmp[, 7], probs = c(0.025, 0.5, 0.975))
  a8[count, ] <- quantile(tmp[, 8], probs = c(0.025, 0.5, 0.975))
  a9[count, ] <- quantile(tmp[, 9], probs = c(0.025, 0.5, 0.975))
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

library(ggplot2)
# Plotting with ggplot2
p1 <- ggplot() +
  geom_point(data = coeff_entropic, aes(x = x, y = median), color = "blue", alpha = 0.5, size = 3) +
  geom_errorbar(data = coeff_entropic, aes(ymin = low_ci, ymax = high_ci, y = median, x = x), width = 0.05, color = "blue", linewidth = 1, alpha = 0.25) +
  geom_point(data = coeff_bme, aes(x = x, y = median), color = "red", alpha = 0.5, size = 3) +
  geom_errorbar(data = coeff_bme, aes(ymin = low_ci, ymax = high_ci, y = median, x = x), width = 0.05, color = "red", linewidth = 1, alpha = 0.25) +
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
  geom_errorbar(data = rf_like, aes(ymin = low_ci, ymax = high_ci, y = median, x = x), width = 0.075, color = "black", size = 1, alpha = 0.25) +
  geom_point(data = rf_entropic, aes(x = x, y = median), color = "blue", alpha = 0.5, size = 3) +
  geom_errorbar(data = rf_entropic, aes(ymin = low_ci, ymax = high_ci, y = median, x = x), width = 0.05, color = "blue", size = 1, alpha = 0.25) +
  geom_point(data = rf_bme, aes(x = x, y = median), color = "red", alpha = 0.5, size = 3) +
  geom_errorbar(data = rf_bme, aes(ymin = low_ci, ymax = high_ci, y = median, x = x), width = 0.05, color = "red", size = 1, alpha = 0.25) +
  geom_vline(aes(xintercept = 0.1), linetype = "dotted", alpha = 0.75) +
  geom_vline(aes(xintercept = 0.5), linetype = "dotted", alpha = 0.75) +
  labs(x = "Rate", y = "1 - Normalised RF distance", title = "Topological accuracy") +
  theme_bw() +
  scale_x_log10() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12), # Adjust the size for tick labels
    plot.title = element_text(size = 16)
  )

p3 <- ggplot() +
  geom_point(data = mae_entropic, aes(x = x, y = 100 * median), color = "blue", alpha = 0.5, size = 3) +
  geom_errorbar(data = mae_entropic, aes(ymin = 100 * low_ci, ymax = 100 * high_ci, y = median, x = x), width = 0.05, color = "blue", size = 1, alpha = 0.25) +
  geom_point(data = mae_bme, aes(x = x, y = 100 * median), color = "red", alpha = 0.5, size = 3) +
  geom_errorbar(data = mae_bme, aes(ymin = 100 * low_ci, ymax = 100 * high_ci, y = median, x = x), width = 0.05, color = "red", size = 1, alpha = 0.25) +
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

library(gridExtra)

grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
