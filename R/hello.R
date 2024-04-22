# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

streakTests <- function(data, succ = 1, fail = 0, n = 1000, conf = .95, dp = 3, pdp = 3,
                        runsT = T, acT = T, lag = 1, cpT = T, prior = 3, graph = T) {
  #
  # Performs tests for streaks in binary data using empirical rather theoretical sampling distributions.
  # Tests are runs and autocorrelation (conditional probabilities to be added later).
  # Accepts data in the form of a vector (single sequence) or a matrix (one sequence per row).
  # For vectors, missing events removed; for matrices, entire row removed if any missing events.
  #
  # Arguments:
  #    data: binary data sequence (or a matrix containing rows of binary data sequences)
  #    succ: numerical code for a "success" event (default = 1)
  #    fail: numerical code for a "failure" event (default = 0)
  #       n: number of bootstrap samples used to generate sampling distributions (default = 1000)
  #    conf: confidence level (default = .95)
  #      dp: number of decimal places (default = 3)
  #     pdp: number of decimal places for p values (default = 3)
  #   runsT: whether to perform the runs test (default = T)
  #     acT: whether to perform the autocorrelation test (default = T)
  #     lag: lag for autocorrelations (default = 1)
  #     cpT: whether to perform the conditional probabilities test (default = T)
  #   prior: number of prior events for conditional probabilities (default = 3; max = 3)
  #   graph: whether to provide visualation of results (default = T)

  # check whether data object is vector or matrix
  type <- 1 * is.vector(data) + 2 * (is.matrix(data))

  # code success and failure as 1 and 0; remove any missing data
  if (type == 1) {
    z <- rep(NA, length(data))
    z[data == succ] <- 1
    z[data == fail] <- 0
    if (sum(is.na(z)) > 0) {
      z <- z[!is.na(z)]
    }
  } else {
    z <- matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2])
    dimnames(z) <- dimnames(data)
    z[data == succ] <- 1
    z[data == fail] <- 0
    i <- 1
    while (i <= dim(z)[1]) {
      if (sum(is.na(z[i, ])) > 0) {
        z <- z[-i, ]
      } else {
        i <- i + 1
      }
    }
  }

  # convert vector to single-row matrix
  if (type == 1) {
    z <- t(as.matrix(z))
  }
  nz <- dim(z)[2]
  nr <- dim(z)[1]

  # embed functions
  runCount <- function(z) {
    return(length(rle(z)$lengths))
  }
  autoCorr <- function(z, lag) {
    n <- length(z)
    x <- z[1:(n - lag)]
    y <- z[(1 + lag):n]
    return(cor(x, y))
  }
  condProbs <- function(z, prior = 3) {
    cps <- matrix(NA, nrow = 1, ncol = 2 * prior + 1)
    sm <- "1"
    fm <- "0"
    for (i in 1:prior) {
      sm <- c(paste(c(rep("0", i), "1"), collapse = ""), sm, paste(c(rep("1", i), "1"), collapse = ""))
      fm <- c(paste(c(rep("0", i), "0"), collapse = ""), fm, paste(c(rep("1", i), "0"), collapse = ""))
    }
    for (i in 1:(2 * prior + 1)) {
      j <- abs(i - (prior + 1))
      if (j == 0) {
        zs <- paste0(z[1:nz])
      }
      if (j == 1) {
        zs <- paste0(z[1:(nz - 1)], z[2:nz])
      }
      if (j == 2) {
        zs <- paste0(z[1:(nz - 2)], z[2:(nz - 1)], z[3:nz])
      }
      if (j == 3) {
        zs <- paste0(z[1:(nz - 3)], z[2:(nz - 2)], z[3:(nz - 1)], z[4:nz])
      }
      sum.s <- sum(zs == sm[i])
      sum.f <- sum(zs == fm[i])
      cps[i] <- sum.s / (sum.s + sum.f)
    }
    return(cps)
  }

  # calculate observed values
  if (runsT) {
    runs.obs <- mean(apply(z, 1, runCount))
  }
  if (acT) {
    ac.obs <- mean(apply(z, 1, autoCorr, lag = lag))
  }
  if (cpT) {
    cp.obs <- rbind(apply(z, 1, condProbs, prior = prior))
    cp.obs <- t(as.matrix(apply(cp.obs, 1, mean)))
    lab <- "p(S)"
    for (i in 1:prior) {
      lab <- c(paste(c("p(S | ", rep("F", i), ")"), collapse = ""), lab, paste(c("p(S | ", rep("S", i), ")"), collapse = ""))
    }
    dimnames(cp.obs) <- list("  Observed", lab)
  }

  # generate sampling distributions
  runs.exp <- ac.exp <- rep(NA, n)
  cp.exp <- matrix(NA, nrow = n, ncol = 2 * prior + 1)
  for (i in 1:n) {
    rseq <- z[, sample(1:nz, replace = T)]
    if (runsT) {
      if (is.vector(rseq)) {
        runs.exp[i] <- length(rle(rseq)$lengths)
      } else {
        runs.exp[i] <- mean(apply(rseq, 1, runCount))
      }
    }
    if (acT) {
      if (is.vector(rseq)) {
        ac.exp[i] <- autoCorr(rseq, lag = lag)
      } else {
        ac.exp[i] <- mean(apply(rseq, 1, autoCorr, lag = lag))
      }
    }
    if (cpT) {
      if (is.vector(rseq)) {
        rseq <- t(as.matrix(rseq))
      }
      cp.temp <- rbind(apply(rseq, 1, condProbs, prior = prior))
      cp.exp[i, ] <- t(as.matrix(apply(cp.temp, 1, mean)))
    }
  }

  # provide text output
  z.crit <- qnorm(1 - (1 - conf) / 2)
  if (runsT) {
    cat("Observed number of runs =", round(runs.obs, dp), "\n")
    cat("  Expected number of runs:\n", sep = "")
    m.exp <- mean(runs.exp)
    se.exp <- sd(runs.exp)
    ci.lower <- m.exp - z.crit * se.exp
    ci.upper <- m.exp + z.crit * se.exp
    cat("    M = ", round(m.exp, dp), ", SE = ", round(se.exp, dp),
        ", CI = ", round(ci.lower, dp), " to ", round(ci.upper, dp), "\n", sep = "")
    z.value <- (runs.obs - m.exp) / se.exp
    p.value <- pnorm(z.value, lower.tail = T)
    cat("  Test of null hypothesis (obs >= exp):\n")
    cat("    Normal approximation:  z = ", round(z.value, dp), ", p = ", round(p.value, pdp), "\n\n", sep = "")
  }
  if (acT) {
    cat("Observed autocorrelation =", round(ac.obs, dp), "\n")
    cat("  Expected autocorrelation:\n", sep = "")
    m.exp <- mean(ac.exp)
    se.exp <- sd(ac.exp)
    ci.lower <- m.exp - z.crit * se.exp
    ci.upper <- m.exp + z.crit * se.exp
    cat("    M = ", round(m.exp, dp), ", SE = ", round(se.exp, dp),
        ", CI = ", round(ci.lower, dp), " to ", round(ci.upper, dp), "\n", sep = "")
    z.value <- (ac.obs - m.exp) / se.exp
    p.value <- pnorm(z.value, lower.tail = F)
    cat("  Test of null hypothesis (obs <= exp):\n")
    cat("    Normal approximation:  z = ", round(z.value, dp), ", p = ", round(p.value, pdp), "\n\n", sep = "")
  }
  if (cpT) {
    cat("Conditional probabilities:\n")
    cp.out <- cp.obs
    cp.out <- rbind(cp.out, apply(cp.exp, 2, mean, na.rm = T))
    cp.out <- rbind(cp.out, apply(cp.exp, 2, sd, na.rm = T))
    cp.out <- rbind(cp.out, cp.out[2, ] - z.crit * cp.out[3, ])
    cp.out <- rbind(cp.out, cp.out[2, ] + z.crit * cp.out[3, ])
    dimnames(cp.out)[[1]] <- c("  Observed", "  Exp M", "  Exp SE", "  CI[lower]", "  CI[upper]")
    print(round(cp.out, dp))
    cat("\n")
  }
  cat("Notes:\n")
  cat("  Number of sequences =", nr, "\n")
  cat("  Sequence length =", nz, "\n")
  cat("  Expected values generated using", n, "random sequences for each observed sequence\n")
  cat("  SE estimated as SD of empirical sampling distribution\n")
  cat("  CI constructed using normal approximation\n")
  cat("  Confidence level for CI = ", round(100 * conf), "%\n", sep = "")
  cat("  p values for 1-tailed tests\n")
  cat("  Number of decimal places =", dp, "\n")
  cat("  Number of decimal places for p values =", pdp, "\n")
  cat("  R drops any 0s at the end of a decimal value\n\n")

  # plot graphs
  if (graph) {
    if (runsT) {
      dev.new(height = 4, width = 6)
      par(cex = .75)
      min.max <- c(min(c(runs.exp, runs.obs)), max(c(runs.exp, runs.obs)))
      runs.hist <- hist(runs.exp, breaks = 35, plot = F)
      y.max <- max(runs.hist$counts)
      xs <- c(sort(runs.exp)[.05 * n], min(runs.exp))
      adj <- .5 + .5 * (runs.obs > min.max[1] + (min.max[2] - min.max[1]) * .9) -
        .5 * (runs.obs < min.max[1] + (min.max[2] - min.max[1]) * .1)
      plot(runs.hist, xlim = min.max, ylim = c(0, y.max * 1.25), main = "Runs Test", xlab = "Number of Runs")
      box(which = "plot")
      lines(x = c(runs.obs, runs.obs), y = c(0, y.max * 1.10), lwd = 2)
      text(x = runs.obs, y = y.max * 1.20, paste("Observed Runs =", round(runs.obs, dp)), adj = adj)
      runs.sig <- sort(runs.exp)[.05 * n]
      lines(x = c(xs[1], xs[1]), y = c(0, y.max), lty = 3)
      lines(x = xs, y = rep(y.max, 2), lty = 3)
      lines(x = c(xs[2], xs[2]), y = c(0, y.max), lty = 3)
      text(x = mean(xs), y = y.max * .8, paste("Most\nStreaky\n5%"))
    }
    if (acT) {
      dev.new(height = 4, width = 6)
      par(cex = .75)
      min.max <- c(min(c(ac.exp, ac.obs)), max(c(ac.exp, ac.obs)))
      ac.hist <- hist(ac.exp, breaks = 35, plot = F)
      y.max <- max(ac.hist$counts)
      xs <- c(sort(ac.exp)[.95 * n], max(ac.exp))
      adj <- .5 + .5 * (ac.obs > min.max[1] + (min.max[2] - min.max[1]) * .9) -
        .5 * (ac.obs < min.max[1] + (min.max[2] - min.max[1]) * .1)
      plot(ac.hist, xlim = min.max, ylim = c(0, y.max * 1.25), main = "Autocorrelations", xlab = "Autocorrelation")
      box(which = "plot")
      lines(x = c(ac.obs, ac.obs), y = c(0, y.max * 1.10), lwd = 2)
      text(x = ac.obs, y = y.max * 1.20, paste("Observed r =", round(ac.obs, dp)), adj = adj)
      lines(x = c(xs[1], xs[1]), y = c(0, y.max), lty = 3)
      lines(x = xs, y = rep(y.max, 2), lty = 3)
      lines(x = c(xs[2], xs[2]), y = c(0, y.max), lty = 3)
      text(x = mean(xs), y = y.max * .8, paste("Most\nStreaky\n5%"))
    }
    if (cpT) {
      dev.new(height = 4, width = 6)
      par(cex = .75)
      xs <- 1:(2 * prior + 1)
      plot(x = xs, y = cp.out[1, ], ylim = c(0, 1), type = "b", xaxt = "n",
           main = "Conditional Probabilities", xlab = "Condition", ylab = "Probability")
      axis(side = 1, at = xs, labels = dimnames(cp.out)[[2]])
      polygon(x = c(xs, rev(xs)), y = c(cp.out[4, xs], cp.out[5, rev(xs)]), col = "gray90", border = NA)
      points(xs, y = cp.out[1, ], pch = 19)
      lines(xs, y = cp.out[1, ], lwd = 2)
      legend(x = mean(xs) + .5, y = .10, legend = c("Observed", "Expected"), bty = "n", lty = c(1, 0),
             lwd = c(2, 0), pch = c(19, 15), pt.cex = c(1, 3), col = c(1, "gray90"), xjust = .5, yjust = .5)
    }
  }
}


# test data
#   x0 = vector of 100 events with no more streaks than expected
#   x1 = vector of 100 events with more streaks than expected (small effect)
#   x2 = vector of 100 events with more streaks than expected (large effect)
#   y0 = matrix of 50 rows, each with 100 events with no more streaks than expected
#   y1 = matrix of 50 rows, each with 100 events with more streaks than expected (small effect)
#   y2 = matrix of 50 rows, each with 100 events with more streaks than expected (large effect)

set.seed(1)
x0 <- rbinom(100, 1, .5)
x1 <- c(rbinom(25, 1, .70), rbinom(25, 1, .30), rbinom(25, 1, .70), rbinom(25, 1, .30))
x2 <- c(rbinom(25, 1, .85), rbinom(25, 1, .15), rbinom(25, 1, .85), rbinom(25, 1, .15))
y0 <- y1 <- y2 <- matrix(0, nrow = 50, ncol = 100)
for (i in 1:50) {
  y0[i,] <- rbinom(100, 1, .5)
  y1[i,] <- c(rbinom(25, 1, .70), rbinom(25, 1, .30), rbinom(25, 1, .70), rbinom(25, 1, .30))
  y2[i,] <- c(rbinom(25, 1, .85), rbinom(25, 1, .15), rbinom(25, 1, .85), rbinom(25, 1, .15))
}
