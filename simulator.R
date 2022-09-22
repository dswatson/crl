### Simulations for contrastive ATE learning ###

# Load libraries, register cores
library(data.table)
library(pcalg)
library(RBGL)
library(doMC)
library(doRNG)
registerDoMC(8)

# Set seed
set.seed(123)

################################################################################

### SIMULATION ###

#' @param n Sample size for the observational regime.
#' @param d Dimensionality.
#' @param n_fctr Sample size for each interventional regime is \code{1/n_fctr}
#'   as large as \code{n}.
#' @param r2 Proportion of variance explained for all structural relationships. 
#' @param lin_pr Probability that an edge denotes a linear relationship.
#' @param conf_pr Probability of confounding between any two edges.
#' @param sp Average sparsity of the graph. Note that this must be high for 
#'   \code{method = "barabasi"} or else you'll run into errors.
#' @param method Method used for generating the graph structure. Options are
#'   \code{"er"} for Erdós-Rényi and \code{"barabasi"} for Barabási-Albert.
#' @param pref Strength of preferential attachment if \code{method = "barabasi"}.
#' 

# Data simulation function
# Note: lower triangular adj_mat means that column is a parent of row
sim_dat <- function(n, d, n_fctr, r2, lin_pr, conf_pr, sp, method, pref) {
  # Optionally apply nonlinear transformations
  prep <- function(dat, pr) {
    out <- dat
    if (pr < 1) {
      # Pick features to transform
      n_nl <- round((1 - pr) * ncol(dat))
      if (n_nl > 0) {
        tmp <- data.table(idx = sample.int(ncol(dat), size = n_nl))
        tmp[, nl := sample(c('sq', 'sqrt', 'sftpls', 'relu'), 
                           size = n_nl, replace = TRUE)]
        out[, tmp[nl == 'sq', idx]] <- dat[, tmp[nl == 'sq', idx]]^2
        out[, tmp[nl == 'sqrt', idx]] <- sqrt(abs(dat[, tmp[nl == 'sqrt', idx]]))
        out[, tmp[nl == 'sftpls', idx]] <- log(1 + exp(dat[, tmp[nl == 'sftpls', idx]]))
        out[, tmp[nl == 'relu', idx]] <- ifelse(dat[, tmp[nl == 'relu', idx]] > 0, 
                                                dat[, tmp[nl == 'relu', idx]], 0)
      }
    }
    return(out)
  }
  # Generate noise
  sim_noise <- function(signal, r2, n) {
    var_mu <- var(signal)
    if (var_mu == 0) var_mu <- 1
    var_noise <- (var_mu - r2 * var_mu) / r2
    noise <- rnorm(n, sd = sqrt(var_noise))
    return(noise)
  }
  # Simulate graph
  m <- (1 - sp) * (d - 1)
  gr <- randDAG(d, m, method = method, par1 = pref, weighted = FALSE)
  t_srt <- as.numeric(tsort(gr))
  x_idx <- data.table(x = seq_len(d), g = t_srt[seq_len(d)])
  # Confounding variables
  x_labs <- paste0('x', seq_len(d))
  conf_idx <- as.data.table(do.call(rbind, combn(x_labs, 2, simplify = FALSE)))
  conf_idx[, conf := rbinom(nrow(conf_idx), size = 1, prob = conf_pr)]
  conf_idx <- conf_idx[conf == 1]
  d_u <- conf_idx[, sum(conf)]
  # Record adjacency matrix 
  if (d_u > 0) {
    u_labs <- paste0('u', seq_len(d_u))
  } else {
    u_labs <- NULL
  }
  adj_mat <- matrix(0, nrow = d_u + d, ncol = d_u + d, 
                    dimnames = list(c(u_labs, x_labs), c(u_labs, x_labs)))
  diag(adj_mat) <- NA_real_
  pa <- lapply(seq_len(d), function(j) {
    which(sapply(seq_len(d), function(i) {
      t_srt[j] %in% gr@edgeL[[i]]$edges
    }))
  })
  for (j in seq_len(d)) {
    p_idx <- x_idx[g %in% pa[[j]], x]
    adj_mat[j + d_u, p_idx + d_u] <- 1
  }
  for (j in seq_len(d_u)) {
    c_idx <- as.character(c(conf_idx[j, V1], conf_idx[j, V2]))
    adj_mat[c_idx, j] <- 1
  }
  # Sample structural weights from Rademacher distribution
  beta <- lapply(seq_len(d + d_u), function(j) {
    sample(c(1, -1), size = sum(adj_mat[j, ], na.rm = TRUE), replace = TRUE)
  })
  # Compute X recursively
  dag_loop <- function(sigma) {
    n_tmp <- ifelse(sigma == 0, n, round(n / n_fctr))
    u <- matrix(rnorm(n_tmp * d_u), ncol = d_u, dimnames = list(NULL, u_labs))
    x <- matrix(nrow = n_tmp, ncol = d, dimnames = list(NULL, x_labs))
    for (j in seq_len(d)) {
      if (sigma == j) {
        x[, j] <- 0
      } else if (sum(adj_mat[j + d_u, ], na.rm = TRUE) == 0) {
        x[, j] <- rnorm(n_tmp)
      } else {
        pa_idx <- rownames(adj_mat)[which(adj_mat[j + d_u, ] == 1)]
        if (any(grepl('u', pa_idx))) {
          pa_u <- as.matrix(u[, colnames(u) %in% pa_idx])
        } else {
          pa_u <- NULL
        }
        if (any(grepl('x', pa_idx))) {
          pa_x <- as.matrix(x[, colnames(x) %in% pa_idx])
        } else {
          pa_x <- NULL
        }
        pa_j <- as.matrix(prep(cbind(pa_u, pa_x), lin_pr))
        signal_x <- as.numeric(pa_j %*% beta[[j + d_u]])
        x[, j] <- signal_x + sim_noise(signal_x, r2, n_tmp)
      }
    }
    out <- data.table(x)[, 'sigma' := sigma]
    return(out)
  }
  out <- foreach(s = 0:d, .combine = rbind) %dorng% dag_loop(s)
  # Export
  params <- list(
    'n' = n, 'd' = d, 'n_fctr' = n_fctr, 
    'r2' = r2, 'lin_pr' = lin_pr, 'conf_pr' = conf_pr, 
    'sp' = sp, 'method' = method, 'pref' = pref
  )
  out <- list('dat' = out, 'adj_mat' = adj_mat, 'beta' = beta, 
              'params' = params)
  return(out)
}















