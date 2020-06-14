run_sim = function(C, percent_vax, strategy, u = u_constant, v_e = v_e_constant, frac_age = age_demo, npop = pop_total, N = N_i, nage = num_groups){

  # Disease Tranmission
  nu    <- 1/3 # incubation period (E -> I), ref: Davies
  gamma <- 1/5 # recovery period (I -> R), ref: Davies

  # _____________________________________________________________________
  # VACCINE STRATEGIES ----
  #     propall: distribute proportionally to each age group
  #     propkids: distribute proportionally to age groups < 20
  #     propadults: distribute proportionally to age groups 20-49
  #     propelderly: distribute proportionally to age groups 60+
  # _____________________________________________________________________
  nvax <- percent_vax*npop 
  
  # Initialize simulation with 1 infected in each age group 
  I_0    <- rep(1,nage)
  E_0    <- rep(0,nage)
  R_0    <- rep(0,nage)
  
  if (strategy == "no vax"){
    S_0    <- N-I_0
    V_0    <- rep(0,nage)
    
    inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, V=V_0)
  } else if (strategy == "all"){
    # prop strategy
    vax_propall <- nvax*frac_age
    vax_propall[vax_propall > N] <- N[vax_propall > N] # don't exceed num people in pop
    
    S_0    <- N-I_0-vax_propall
    V_0    <- vax_propall
    
    inits <- c(S=S_0, E=E_0, I=I_0, R=R_0, V=V_0)
  } else if (strategy == "kids"){
    # propkids strategy
    nkids <- N[1] + N[2]
    vax_dist_propkids <- rep(0, nage)
    vax_dist_propkids[1] <- N[1]/nkids
    vax_dist_propkids[2] <- N[2]/nkids
    
    vax_propkids <- nvax*vax_dist_propkids
    vax_propkids[vax_propkids > N] <- N[vax_propkids > N]
    
    S_0    <- N-I_0-vax_propkids
    V_0    <- vax_propkids
    
    inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, V=V_0)
  } else if (strategy == "adults") {
    # propadults strategy
    nadults <- N[3] + N[4] + N[5]
    vax_dist_propadults <- rep(0, nage)
    vax_dist_propadults[3] <- N[3]/nadults
    vax_dist_propadults[4] <- N[4]/nadults
    vax_dist_propadults[5] <- N[5]/nadults
    
    vax_propadults <- nvax*vax_dist_propadults
    vax_propadults[vax_propadults > N] <- N[vax_propadults > N]
    
    S_0    <- N-I_0-vax_propadults
    V_0    <- vax_propadults
    
    inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, V=V_0)
  } else if (strategy == "elderly") {
    # propelderly strategy
    nelderly <- N[7] + N[8] + N[9]
    vax_dist_propelderly <- rep(0, nage)
    vax_dist_propelderly[7] <- N[7]/nelderly
    vax_dist_propelderly[8] <- N[8]/nelderly
    vax_dist_propelderly[9] <- N[9]/nelderly
    
    vax_propelderly <- nvax*vax_dist_propelderly
    vax_propelderly[vax_propelderly > N] <- N[vax_propelderly > N]
    
    S_0    <- N-I_0-vax_propelderly
    V_0    <- vax_propelderly
    
    inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, V=V_0)
  }
  # _____________________________________________________________________
  # NUMERICALLY SOLVE ----
  # _____________________________________________________________________
  parameters = list(beta=u, nu=nu, gamma=gamma, C=C, v_e=v_e)
  
  t <- seq(0,150,1) 
  
  df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
  df
}
