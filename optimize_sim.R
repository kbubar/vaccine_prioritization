optimize_sim = function(vax_vec) {
  # INPUT:
  #    vax_vector: vector of initial vaccine allocation
  
  # parameters
  C       <- C                      
  u       <- u_constant
  v_e     <- v_e_constant
  
  # constants
  frac_age<- age_demo
  npop    <- pop_total
  N       <- N_i
  nage    <- num_groups
  nu      <- 1/3 # incubation period (E -> I), ref: Davies
  gamma   <- 1/5 # recovery period (I -> R), ref: Davies
  
  # Initialize simulation with 1 infected in each age group 
  I_0    <- rep(1,nage)
  E_0    <- rep(0,nage)
  R_0    <- rep(0,nage)
  
  S_0    <- N-I_0-vax_vec
  V_0    <- vax_vec
    
  inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, V=V_0)

  # _____________________________________________________________________
  # NUMERICALLY SOLVE ----
  # _____________________________________________________________________
  parameters = list(beta=u, nu=nu, gamma=gamma, C=C, v_e=v_e)
  
  t <- seq(0,150,1) 
  
  df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))

  if (to_minimize == "cases"){
    tot <- compute_total_cases(df)
  } else if (to_minimize == "deaths"){
    tot <- compute_total_deaths(df)
  }
  
  # return(df)
  return(tot)
}
