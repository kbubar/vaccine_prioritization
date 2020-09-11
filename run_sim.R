run_sim = function(C, percent_vax, strategy, u = u_constant, v_e = v_e_constant, 
                   frac_age = age_demo, npop = pop_total, N = N_i, nage = num_groups, 
                   sero = sero_none, sero_testing = FALSE){

  # Disease Tranmission
  d_E    <- 1/3 # incubation period (E -> I), ref: Davies
  d_I <- 1/5 # recovery period (I -> R), ref: Davies

  # _____________________________________________________________________
  # INITIAL CONDITIONS ----
  # Vaccine strategies distribute proportionally to 
  #     all: all age groups
  #     kids: age groups < 20
  #     adults: age groups 20-49
  #     elderly: age groups 60+
  #     20+: age groups 20+
  # _____________________________________________________________________
  E_0    <- rep(0,nage)
  R_0    <- N * sero
  
  # specify group to vaccinate according to allocation strategy
  if (strategy == "no vax"){
    V_0 <- rep(0, nage)
  } else {
      if (strategy == "all"){ 
      groups <- 1:9
    } else if (strategy == "kids"){ 
      groups <- 1:2
    } else if (strategy == "adults") { 
      groups <- 3:5
    } else if (strategy == "elderly") {
      groups <- 7:9
    } else if (strategy == "20+") {
      groups <- 3:9
    }
    people_to_vax <- sum(N[groups])
    
    vax_proportion <- rep(0, nage)
    vax_proportion[groups] <- N[groups]/people_to_vax
    
    nvax <- percent_vax*npop 
    vax_distribution <- nvax*vax_proportion
    vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
    
    prob_vaccinated <- vax_distribution/N
    
    if (sero_testing == FALSE){
      # account for people who were vaccinated and seropositive
      V_0    <- vax_distribution - N*prob_vaccinated*sero
      #V_0 <- vax_distribution*v_e # all-or-nothing vaccine
    } else if (sero_testing == TRUE){
      # assumes everyone vaccinated was seronegative
      V_0    <- vax_distribution
      temp <- N - R_0
      V_0[V_0 > temp] <- temp[V_0 > temp]
    }
  }
  # initial I: 1 in each age group unless everyone is vaccinated and/or sero positive
  I_0 <- rep(0, nage)
  I_0[(N - R_0 - V_0) > 1] <- 1
  
  S_0 <- N - I_0 - V_0 - R_0
  
  inits <- c(S_0,E_0,I_0,R_0,V_0)
  
  # _____________________________________________________________________
  # NUMERICALLY SOLVE ----
  # _____________________________________________________________________
  parameters = list(u=u, d_E=d_E, d_I=d_I, C=C, v_e=v_e)

  # t <- seq(0,800,1)
  # tfinal <- 800
  # 
  # df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))

  running = TRUE
  t <- c(0:50)
  df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
  t <- t + 50

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
      t <- t + 50
      }
  }

  names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
                    "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
                    "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
                    "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
                    "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")

  return(df)
}

run_sim_nontransmissionblocking = function(C, percent_vax, strategy, alpha, omega, u = u_constant, v_e = v_e_constant, frac_age = age_demo, npop = pop_total, N = N_i, nage = num_groups, sero = sero_none, sero_testing = FALSE){
  
  # Disease Tranmission
  d_E    <- 1/3 # incubation period (E -> I), ref: Davies
  d_I <- 1/5 # recovery period (I -> R), ref: Davies
  
  # _____________________________________________________________________
  # INITIAL CONDITIONS ----
  # Vaccine strategies distribute proportionally to 
  #     all: all age groups
  #     kids: age groups < 20
  #     adults: age groups 20-49
  #     elderly: age groups 60+
  #     20+: age groups 20+
  # _____________________________________________________________________
  E_0    <- rep(0,nage)
  R_0    <- N * sero
  
  # specify group to vaccinate according to allocation strategy
  if (strategy == "no vax"){
    V_0 <- rep(0, nage)
  } else {
    if (strategy == "all"){ 
      groups <- 1:9
    } else if (strategy == "kids"){ 
      groups <- 1:2
    } else if (strategy == "adults") { 
      groups <- 3:5
    } else if (strategy == "elderly") {
      groups <- 7:9
    } else if (strategy == "20+") {
      groups <- 3:9
    }
    people_to_vax <- sum(N[groups])
    
    vax_proportion <- rep(0, nage)
    vax_proportion[groups] <- N[groups]/people_to_vax
    
    nvax <- percent_vax*npop 
    vax_distribution <- nvax*vax_proportion
    vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
    
    prob_vaccinated <- vax_distribution/N
    
    if (sero_testing == FALSE){
      # account for people who were vaccinated and seropositive
      #V_0    <- vax_distribution - N*prob_vaccinated*sero
      V_0 <- vax_distribution*v_e # all or nothing vaccine
    } else if (sero_testing == TRUE){
      # assumes everyone vaccinated was seronegative
      V_0    <- vax_distribution
      temp <- N - R_0
      V_0[V_0 > temp] <- temp[V_0 > temp]
    }
  }
  # initial I: 1 in each age group unless everyone is vaccinated and/or sero positive
  I_0 <- rep(0, nage)
  I_0[(N - R_0 - V_0) > 1] <- 1
  
  S_0 <- N - I_0 - V_0 - R_0
  
  inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, SV=V_0, EV = rep(0,nage), IV = rep(0,nage), RV = rep(0,nage))
  
  # _____________________________________________________________________
  # NUMERICALLY SOLVE ----
  # _____________________________________________________________________
  parameters = list(u=u, d_E = d_E, d_I = d_I, C=C, v_e=v_e, alpha = alpha, omega = omega)

  t <- seq(0,1000,1)
  tfinal <- 1000
  
  df <- as.data.frame(lsoda(inits, t, calculate_derivatives_nontransmissionblocking, parameters))
  
  I_tfinal <- sum(df[(tfinal+1),20:28])
  
  if (I_tfinal > 1){
    print("Warning: simulation may not have run long enough")
  }
  
  df
}

run_sim_later_vax <- function(t_infected, strategy, percent_vax, sero_testing = FALSE, 
                              npop = pop_total, N = N_i, nage = num_groups) {
  # t_infected = time to distribute vaccine wrt % of infections that have occured
  
  df_baseline <- list_all[[1]]
  R_time <- apply(df_baseline[29:37], 1, sum)
  R_time <- (R_time/pop_total)*100
  tot_infected <- compute_total_cases(list_all[[1]])
  t_init <- length(R_time[R_time <= tot_infected*t_infected])
  
  inits <- as.numeric(df_baseline[t_init, -(1)])
  S_0 <- inits[1:9]
  E_0 <- inits[10:18]
  I_0 <- inits[19:27]
  R_0 <- inits[28:36]
  
  sero <- R_0/N_i
  # age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  # df <- data.frame(age_groups, sero)
  # 
  # theme_set(theme_minimal(base_size = 15))
  # 
  #  ggplot(df, aes(y=sero*100, x=age_groups)) +
  #    geom_bar(position="stack", stat="identity") +
  #    xlab("Age group") +
  #    ylab("Percent Seropositive") +
  #    scale_x_discrete(breaks=c("0-9","20-29", "40-49", "60-69", "80+")) +
  #    theme(plot.subtitle = element_text(size = 9, face = "italic", hjust = 1),
  #          plot.title = element_text(hjust = 0.5)) +
  #    ylim(0, 65)
  
          
  if (strategy == "all"){ 
    groups <- 1:9
  } else if (strategy == "kids"){ 
    groups <- 1:2
  } else if (strategy == "adults") { 
    groups <- 3:5
  } else if (strategy == "elderly") {
    groups <- 7:9
  } else if (strategy == "20+") {
    groups <- 3:9
  }
  people_to_vax <- sum(N[groups])
  
  vax_proportion <- rep(0, nage)
  vax_proportion[groups] <- N[groups]/people_to_vax
  
  nvax <- percent_vax*npop 
  vax_distribution <- nvax*vax_proportion
  vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
  
  prob_vaccinated <- vax_distribution/N
  
  if (sero_testing == FALSE){
    # account for people who were vaccinated and seropositive
    V_0    <- vax_distribution - N*prob_vaccinated*sero
    V_0[V_0 > S_0] <- S_0[V_0 > S_0]
  } else if (sero_testing == TRUE){
    # assumes everyone vaccinated was seronegative
    V_0    <- vax_distribution
    V_0[V_0 > S_0] <- S_0[V_0 > S_0]
  }
  S_0 <- S_0 - V_0
  
  inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, V=V_0)
  
  nu    <- 1/3 # incubation period (E -> I), ref: Davies
  gamma <- 1/5 # recovery period (I -> R), ref: Davies
  parameters = list(beta=u_constant, nu=nu, gamma=gamma, C=C, v_e=v_e_constant)

  t <- seq(0,150,1)

  df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
  df
}
