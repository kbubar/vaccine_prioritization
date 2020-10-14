optimize_sim = function(vax_vec) {
  # INPUT:
  #    vax_vector: vector of initial vaccine allocation
  
  # parameters
  C       <- C                    
  u       <- u_var
  v_e     <- v_e_constant
  
  # constants
  frac_age<- age_demo
  npop    <- pop_total
  N       <- N_i
  sero    <- sero_none # sero_belgium
  nage    <- num_groups
  d_E      <- 1/3 # incubation period (E -> I), ref: Davies
  d_I   <- 1/5 # recovery period (I -> R), ref: Davies
  
  E_0    <- rep(0,nage)
  #R_0    <- rep(0,nage)
  R_0    <- N * sero
  V_0    <- vax_vec
  
  # initial I: 1 in each age group unless everyone is vaccinated and/or sero positive
  I_0 <- rep(0, nage)
  I_0[(N - R_0 - V_0) > 1] <- 1
  
  S_0    <- N-I_0-V_0
 
  inits <- c(S_0,E_0,I_0,R_0,V_0)

  # _____________________________________________________________________
  # NUMERICALLY SOLVE ----
  # _____________________________________________________________________
  parameters = list(u=u, d_E=d_E, d_I=d_I, C=C, v_e=v_e)
  
  running = TRUE
  t <- c(0:20)
  df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
  t <- t + 20
  
  while(running == TRUE){
    inits <- as.numeric(df[t[1]+1, -(1)])
    temp <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
    row.names(temp) <- t+1
    temp <- temp[-(1),]
    
    df <- rbind(df, temp)
    
    I_tfinal <- sum(df[(t+1),20:28])
    if (I_tfinal < 1){
      running = FALSE
    } else {
      t <- t + 20
    }
  }
  
  names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
                 "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
                 "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
                 "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
                 "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
  
  if (to_minimize == "cases"){
    tot <- compute_total_cases(df)
  } else if (to_minimize == "deaths"){
    tot <- compute_total_deaths(df)
  }
  
  # return(df)
  return(tot)
}
