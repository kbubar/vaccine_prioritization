setwd("~/Vaccine Strategy/Vaccine_Allocation_Project")

df <- read.csv("BEL_prev.csv", header = FALSE)

N <- dim(df)[1] # num of MCMC trials

prev <- matrix(nrow = 100, ncol = 9)

# Sample from empirical CDF
for (i in 1:9){
  U <- runif(100, 0, 1)
  I <- ceiling(N*U)
  prev_vec <- df[I, i]
  prev[, i] <- prev_vec
}

saveRDS(prev, "BEL_prev_samples.RData")
