# Vectorized simulation of state probabilities ---------------------------------
mlogit_probs_vec <- function(x, beta) {
  n_samples <- nrow(beta[[1]])
  prob <- matrix(NA, 
                 ncol = length(beta) + 1,
                 nrow = nrow(x) * n_samples)
  
  # Reference category
  exb <- sapply(beta, function(beta) exp(x %*% t(beta)))
  exb_sum <- rowSums(exb)
  prob[, 1] <- 1/(1 + exb_sum)
  
  # Remaining categories
  for (i in 1:length(beta)) {
    prob[, i + 1] <- prob[1] * exp(x %*% t(beta[[i]]))
  }
  prob
}

tp_sick_vec <- function(input_data, params) {
  input_data <- data.frame(input_data) # Since using data.table with hesim
  beta <- params[c("sick_sicker", "sick_death")]
  x <- make_x(input_data, beta[[1]])
  mlogit_probs_vec(x, beta)
}

tp_sicker_vec <- function(input_data, params) {
  input_data <- data.frame(input_data) # Since using data.table with hesim
  beta <- params[["sicker_death"]]
  x <- make_x(input_data, beta)
  rate <- c(exp(x %*% t(beta)))
  prob_death <- 1 - exp(-rate)
  as.matrix(
    data.frame(sick = 0,
               sicker = 1 - prob_death,
               death = prob_death)
  )
}

tpmatrix_vec <- function(input_data, params) {
  # Transition probabilities
  tp_sick <- tp_sick_vec(input_data, params)
  tp_sicker <- tp_sicker_vec(input_data, params)
  tpmat_vec <- hesim::tpmatrix(
    tp_sick,
    tp_sicker,
    0, 0, 1
  )
  colnames(tpmat_vec) <- hesim::tpmatrix_names(
    states = c("sick", "sicker", "death"), prefix = "", sep = "."
  )
  
  # Row IDs
  n_samples <- nrow(params[[1]])
  tpmat_id <- hesim::tpmatrix_id(input_data, n_samples = n_samples)
  
  # Return
  list(tpmat = tpmat_vec,
       id = tpmat_id)
}