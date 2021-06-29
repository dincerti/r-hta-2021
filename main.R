rm(list = ls())
set.seed(10)
N_SAMPLES <- 100
N_CYCLES <- 20

# Non vectorized approach ------------------------------------------------------
source("functions.R")
input_data <- make_input_data()
params <- set_params(input_data)
sim_stateprobs1(n_cycles = N_SAMPLES, input_data[1, ], params)
system.time(
  out <- sim_stateprobs(input_data, params, n_cycles = 20, 
                        n_samples = N_SAMPLES)
)

# Vectorized approach with hesim -----------------------------------------------
library("hesim")
source("functions-vectorized.R")
input_data2 <- expand(
  hesim::hesim_data(
    strategies = make_strategies(),
    patient = make_patients()
  )
)

# Run/time vectorized version
ptm <- proc.time()
params_sample <- sample_params(n = N_SAMPLES, params)
tpmat_vec <- tpmatrix_vec(input_data2, params_sample)
trans_model <- hesim::CohortDtstmTrans$new(
  params = hesim::tparams_transprobs(tpmat_vec$tpmat,
                                     tpmat_vec$id)
)
out2 <- trans_model$sim_stateprobs(n_cycles = N_CYCLES)
proc.time() - ptm

# Output to check --------------------------------------------------------------
tpmatrix(input_data[1, ], sample_params(n = 1, params))

# Output to share on slides ----------------------------------------------------
input_data[c(1, 2, 49, 50), ]
params
params_sample$sick_sicker[1, ]
params_sample$sicker_death[1, ]
out[1:22, ]
out2