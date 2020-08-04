optimize_sim = function(vax_vec) {
  # INPUT:
  #    vax_vector: vector of initial vaccine allocation
  
  # parameters
  C       <- C_belgium                     
  u       <- u_var
  v_e     <- v_e_constant
  
  # constants
  frac_age<- age_demo_belgium
  npop    <- pop_total
  N       <- N_i
  sero    <- sero_belgium
  nage    <- num_groups
  nu      <- 1/3 # incubation period (E -> I), ref: Davies
  gamma   <- 1/5 # recovery period (I -> R), ref: Davies
  
  E_0    <- rep(0,nage)
  #R_0    <- rep(0,nage)
  R_0    <- N * sero
  V_0    <- vax_vec
  
  # initial I: 1 in each age group unless everyone is vaccinated and/or sero positive
  I_0 <- rep(0, nage)
  I_0[(N - R_0 - V_0) > 1] <- 1
  
  S_0    <- N-I_0-V_0
 
  inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, V=V_0)

  # _____________________________________________________________________
  # NUMERICALLY SOLVE ----
  # _____________________________________________________________________
  parameters = list(beta=u, nu=nu, gamma=gamma, C=C, v_e=v_e)
  
  t <- seq(0,1800,1) 
  
  df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
  
  I_tfinal <- sum(df[(1800+1),20:28])
  
  if (I_tfinal > 1){
    print("Warning: simulation may not have run long enough")
  }
  
  if (to_minimize == "cases"){
    tot <- compute_total_cases(df)
  } else if (to_minimize == "deaths"){
    tot <- compute_total_deaths(df)
  }
  
  # return(df)
  return(tot)
}
