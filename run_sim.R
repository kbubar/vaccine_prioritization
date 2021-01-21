run_sim = function(C, percent_vax, strategy, num_perday, v_e_type, v_e = v_e_constant,
                       u = u_var, sero = sero_none, sp = 1, se = 0, refuse_vax = 0.3, syn_sero_compartments = NA){
  # IC: 0.25% of each age group exposted and infected if numperday != 1
  #     Sero+ start in recovered
  #     refuse_vax % proportion of each compartment start in *x (Sx, Ex, Ix and Rx)
  #     everyone else starts in S
  # Runtime: 365 days post start of rollout
  # Rollout: Continuous at num_perday% of total pop/day until all vaccines are distributed
  #          so num_perday = 1 implies anticipatory rollout
  
  E_0 <- Ev_0 <- Ex_0 <- Sv_0 <- Sx_0 <- Iv_0 <- Rv_0 <- D_0 <- rep(0, num_groups)
  R_0 <- (N_i * sero)*(1-refuse_vax)
  Rx_0 <- (N_i * sero)*refuse_vax
  
  if (num_perday == 1){
     I_0 <- rep(1, num_groups)
     Ix_0 <- rep(1, num_groups)
  } else {
     I_0 <- (N_i*0.0025)*(1-refuse_vax) # 0.5% of each age group starts in exposed and infected
     Ix_0 <- (N_i*0.0025)*refuse_vax
     E_0 <- (N_i*0.0025)*(1-refuse_vax)
     Ex_0 <- (N_i*0.0025)*refuse_vax
  }
  
  S_0 <- (N_i - I_0 - E_0 - R_0 - Ix_0 - Ex_0 - Rx_0)*(1-refuse_vax)
  Sx_0 <- (N_i - I_0 - E_0 - R_0 - Ix_0 - Ex_0 - Rx_0)*refuse_vax
  
  # specify group to vaccinate according to allocation strategy
  if (strategy == "all"){ 
    groups <- 1:9
  } else if (strategy == "kids"){ 
    groups <- 1:2
  } else if (strategy == "adults") { 
    groups <- 3:5
  } else if (strategy == "elderly") {
    groups <- 7:9
  } else if (strategy == "twentyplus") {
    groups <- 3:9
  }
  people_to_vax <- sum(S_0[groups] + E_0[groups] + R_0[groups])
  
  vax_proportion <- rep(0, num_groups)
  vax_proportion[groups] <- (S_0[groups] + E_0[groups] + R_0[groups])/people_to_vax
  
  vax_supply <- percent_vax*pop_total
  
  # anticipatory rollout IC - vaccinate people at t=0
  if (num_perday == 1){
    nvax <- vax_supply
    vax_distribution <- nvax*vax_proportion
    S <- S_0
    E <- E_0
    R <- R_0
    vax_eligible <- (S+E)*sp + R*(1-se)
    
    if (any(vax_distribution > vax_eligible)){
      # make sure everyone in the specificed age groups are vaccinated
      if (!all(vax_distribution[groups] > vax_eligible[groups])){
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        
        full_groups <- 1:9
        full_groups <- full_groups[vax_distribution > vax_eligible]
        leftover_groups <- groups[!groups %in% full_groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
        vax_proportion <- rep(0, num_groups)
        vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
        
        vax_leftover_dist <- leftover_vax*vax_proportion
        vax_distribution <- temp + vax_leftover_dist
      }
      #vax_distribution[vax_distribution > (S+E+R)] <- S[vax_distribution > (S+E+R)] + E[vax_distribution > (S+E+R)] + R[vax_distribution > (S+E+R)]
      #distribute vaccines to all other age groups if there's doses left over after strategy specified
      #age groups
      if (any(vax_distribution > vax_eligible)){
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        
        leftover_groups <- 1:9
        leftover_groups <- leftover_groups[!leftover_groups %in% groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
        vax_proportion <- rep(0, num_groups)
        vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
        
        vax_leftover_dist <- leftover_vax*vax_proportion
        vax_distribution <- temp + vax_leftover_dist
      }
    }
    
    alpha <- as.matrix(vax_distribution)/(as.matrix(vax_eligible))
    alpha[alpha == Inf] <- 0
    alpha[is.nan(alpha)] <- 0
    
    if (v_e_type == "leaky") {
      # first move vaccinated people
      S_0  <- S - as.matrix(S*alpha*sp) - as.matrix(S*alpha*(1-sp))
      Sv_0 <- as.matrix(S*alpha*sp)
      Sx_0 <- Sx_0 + as.matrix(S*alpha*(1-sp))
      
    } else {
      S_0  <- S - as.matrix(S*alpha*sp) - as.matrix(S*alpha*(1-sp))
      Sv_0 <- as.matrix(S*alpha*sp*v_e)
      Sx_0 <- Sx_0 + as.matrix(S*alpha*sp*(1-v_e)) + as.matrix(S*alpha*(1-sp)) 
    }
    R_0  <- R - as.matrix(R*alpha*(1-se)) - as.matrix(R*alpha*se)
    Rv_0 <- as.matrix(R*alpha*(1-se))
    Rx_0 <- Rx_0 + as.matrix(R*alpha*se)
    
    vax_supply <- 0
  }
  
  compartments_initial <- c(S_0,Sv_0,Sx_0,E_0,Ev_0,Ex_0,I_0,Iv_0,Ix_0,R_0,Rv_0,Rx_0,D_0,vax_supply)
  
  if (length(syn_sero_compartments) > 1){
    compartments_initial <- c(as.numeric(syn_sero_compartments), vax_supply)
  }
  
  parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_e_type = v_e_type, num_groups = num_groups, 
                    N_i = N_i, num_perday=num_perday, vax_proportion=vax_proportion, groups=groups, sp=sp, se=se, pop_total=pop_total)

  t <- seq(from=0, to=365, by=1) 
  event_times <- seq(from=0, to=floor(vax_supply/(num_perday*pop_total)), by=1)
  df <- as.data.frame(deSolve::lsoda(y=compartments_initial, times=t, func=calculate_derivatives, parms=parameters, events=list(func=move_vaccinated_event, time=event_times)))
  
  names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
                 "Sv1", "Sv2", "Sv3", "Sv4", "Sv5", "Sv6", "Sv7", "Sv8", "Sv9",
                 "Sx1", "Sx2", "Sx3", "Sx4", "Sx5", "Sx6", "Sx7", "Sx8", "Sx9",
                 "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
                 "Ev1", "Ev2", "Ev3", "Ev4", "Ev5", "Ev6", "Ev7", "Ev8", "Ev9",
                 "Ex1", "Ex2", "Ex3", "Ex4", "Ex5", "Ex6", "Ex7", "Ex8", "Ex9",
                 "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
                 "Iv1", "Iv2", "Iv3", "Iv4", "Iv5", "Iv6", "Iv7", "Iv8", "Iv9",
                 "Ix1", "Ix2", "Ix3", "Ix4", "Ix5", "Ix6", "Ix7", "Ix8", "Ix9",
                 "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
                 "Rv1", "Rv2", "Rv3", "Rv4", "Rv5", "Rv6", "Rv7", "Rv8", "Rv9",
                 "Rx1", "Rx2", "Rx3", "Rx4", "Rx5", "Rx6", "Rx7", "Rx8", "Rx9",
                 "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9")
  
  return(df)
}

run_sim_NTB = function(C, percent_vax, strategy, num_perday, ve_S, ve_I, ve_P,
                       u = u_var, sero = sero_none,  sp = 1, se = 0, refuse_vax = 0.3, syn_sero_compartments = NA){
  # New IC: start w/ 0.5% of each age group in I
  # Only run sim for 365 days post start of rollout
  # Vaccine rollout is continuous at 1% of total pop/day until all vaccines are distributed
  
  E_0 <- Ev_0 <- Ex_0 <- Sv_0 <- Sx_0 <- Iv_0 <- Rv_0 <- D_0 <- rep(0, num_groups)
  R_0 <- (N_i * sero)*(1-refuse_vax)
  Rx_0 <- (N_i * sero)*refuse_vax
  
  if (num_perday == 1){
    I_0 <- rep(1, num_groups)
    Ix_0 <- rep(1, num_groups)
  } else {
    I_0 <- (N_i*0.0025)*(1-refuse_vax) # 0.5% of each age group starts in exposed and infected
    Ix_0 <- (N_i*0.0025)*refuse_vax
    E_0 <- (N_i*0.0025)*(1-refuse_vax)
    Ex_0 <- (N_i*0.0025)*refuse_vax
  }
  
  S_0 <- (N_i - I_0 - E_0 - R_0 - Ix_0 - Ex_0 - Rx_0)*(1-refuse_vax)
  Sx_0 <- (N_i - I_0 - E_0 - R_0 - Ix_0 - Ex_0 - Rx_0)*refuse_vax
  
  # specify group to vaccinate according to allocation strategy
  if (strategy == "all"){ 
    groups <- 1:9
  } else if (strategy == "kids"){ 
    groups <- 1:2
  } else if (strategy == "adults") { 
    groups <- 3:5
  } else if (strategy == "elderly") {
    groups <- 7:9
  } else if (strategy == "twentyplus") {
    groups <- 3:9
  }
  people_to_vax <- sum(S_0[groups] + E_0[groups] + R_0[groups])
  
  vax_proportion <- rep(0, num_groups)
  vax_proportion[groups] <- (S_0[groups] + E_0[groups] + R_0[groups])/people_to_vax
  
  vax_supply <- percent_vax*pop_total
  
  # anticipatory rollout IC - vaccinate people at t=0
  if (num_perday == 1){
    nvax <- vax_supply
    vax_distribution <- nvax*vax_proportion
    S <- S_0
    E <- E_0
    R <- R_0
    vax_eligible <- S+E+R
    
    if (any(vax_distribution > vax_eligible)){
      # make sure everyone in the specificed age groups are vaccinated
      if (!all(vax_distribution[groups] > vax_eligible[groups])){
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        
        full_groups <- 1:9
        full_groups <- full_groups[vax_distribution > vax_eligible]
        leftover_groups <- groups[!groups %in% full_groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
        vax_proportion <- rep(0, num_groups)
        vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
        
        vax_leftover_dist <- leftover_vax*vax_proportion
        vax_distribution <- temp + vax_leftover_dist
      }
      #vax_distribution[vax_distribution > (S+E+R)] <- S[vax_distribution > (S+E+R)] + E[vax_distribution > (S+E+R)] + R[vax_distribution > (S+E+R)]
      #distribute vaccines to all other age groups if there's doses left over after strategy specified
      #age groups
      if (any(vax_distribution > vax_eligible)){
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        
        leftover_groups <- 1:9
        leftover_groups <- leftover_groups[!leftover_groups %in% groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
        vax_proportion <- rep(0, num_groups)
        vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
        
        vax_leftover_dist <- leftover_vax*vax_proportion
        vax_distribution <- temp + vax_leftover_dist
      }
    }
    
    alpha <- as.matrix(vax_distribution)/(as.matrix(vax_eligible))
    alpha[alpha == Inf] <- 0
    alpha[is.nan(alpha)] <- 0
    
    # first move vaccinated people
    S_0  <- S - as.matrix(S*alpha)
    Sv_0 <- as.matrix(S*alpha)
      
    R_0  <- R - as.matrix(R*alpha)
    Rv_0 <- as.matrix(R*alpha)
    
    vax_supply <- 0
  }
  
  compartments_initial <- c(S_0,Sv_0,Sx_0,E_0,Ev_0,Ex_0,I_0,Iv_0,Ix_0,R_0,Rv_0,Rx_0,D_0,vax_supply)
  
  if (length(syn_sero_compartments) > 1){
    compartments_initial <- c(as.numeric(syn_sero_compartments), vax_supply)
  }
  
  parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e_type = "leaky", num_groups = num_groups, N_i = N_i, 
                    sp = sp, se = se, num_perday=num_perday, vax_proportion=vax_proportion, groups=groups, pop_total=pop_total,
                    ve_S = ve_S, ve_I = ve_I, ve_P = ve_P)  
  t <- seq(from=0, to=365, by=1) 
  event_times <- seq(from=0, to=floor(vax_supply/(num_perday*pop_total)), by=1)
  df <- as.data.frame(deSolve::lsoda(y=compartments_initial, times=t, func=calculate_derivatives_NTB, parms=parameters, 
                                     events=list(func=move_vaccinated_event, time=event_times)))
  
  names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
                 "Sv1", "Sv2", "Sv3", "Sv4", "Sv5", "Sv6", "Sv7", "Sv8", "Sv9",
                 "Sx1", "Sx2", "Sx3", "Sx4", "Sx5", "Sx6", "Sx7", "Sx8", "Sx9",
                 "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
                 "Ev1", "Ev2", "Ev3", "Ev4", "Ev5", "Ev6", "Ev7", "Ev8", "Ev9",
                 "Ex1", "Ex2", "Ex3", "Ex4", "Ex5", "Ex6", "Ex7", "Ex8", "Ex9",
                 "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
                 "Iv1", "Iv2", "Iv3", "Iv4", "Iv5", "Iv6", "Iv7", "Iv8", "Iv9",
                 "Ix1", "Ix2", "Ix3", "Ix4", "Ix5", "Ix6", "Ix7", "Ix8", "Ix9",
                 "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
                 "Rv1", "Rv2", "Rv3", "Rv4", "Rv5", "Rv6", "Rv7", "Rv8", "Rv9",
                 "Rx1", "Rx2", "Rx3", "Rx4", "Rx5", "Rx6", "Rx7", "Rx8", "Rx9",
                 "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9")
  
  return(df)
}
