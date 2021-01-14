# run_sim = function(C, percent_vax, strategy, u = u_constant, v_e = v_e_constant, 
#                    frac_age = age_demo, npop = pop_total, N = N_i, nage = num_groups, 
#                    sero = sero_none, sero_testing = FALSE){
# 
#   # Disease Tranmission
#   d_E    <- 1/3 # incubation period (E -> I), ref: Davies
#   d_I <- 1/5 # recovery period (I -> R), ref: Davies
# 
#   # _____________________________________________________________________
#   # INITIAL CONDITIONS 
#   # Vaccine strategies distribute proportionally to 
#   #     all: all age groups
#   #     kids: age groups < 20
#   #     adults: age groups 20-49
#   #     elderly: age groups 60+
#   #     20+: age groups 20+
#   # _____________________________________________________________________
#   E_0    <- rep(0,nage)
#   R_0    <- N * sero
#   
#   # specify group to vaccinate according to allocation strategy
#   if (strategy == "no vax"){
#     V_0 <- rep(0, nage)
#   } else {
#       if (strategy == "all"){ 
#       groups <- 1:9
#     } else if (strategy == "kids"){ 
#       groups <- 1:2
#     } else if (strategy == "adults") { 
#       groups <- 3:5
#     } else if (strategy == "elderly") {
#       groups <- 7:9
#     } else if (strategy == "20+") {
#       groups <- 3:9
#     }
#     people_to_vax <- sum(N[groups])
#     
#     vax_proportion <- rep(0, nage)
#     vax_proportion[groups] <- N[groups]/people_to_vax
#     
#     nvax <- percent_vax*npop 
#     vax_distribution <- nvax*vax_proportion
#     vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
#     
#     prob_vaccinated <- vax_distribution/N
#     
#     if (sero_testing == FALSE){
#       # account for people who were vaccinated and seropositive
#       #V_0    <- vax_distribution - N*prob_vaccinated*sero # leaky vaccine 
#       V_0 <- vax_distribution*v_e - N*prob_vaccinated*v_e*sero # all-or-nothing vaccine
#     } else if (sero_testing == TRUE){
#       # assumes everyone vaccinated was seronegative
#       #V_0    <- vax_distribution # leaky vaccine
#       V_0    <- vax_distribution*v_e # all_or_nothing vaccine
#       temp <- N - R_0
#       V_0[V_0 > temp] <- temp[V_0 > temp]
#     }
#     
#     #Extend strategy after everyone in groups has been vaccinated
#     percent_strat <- sum(age_demo[groups]) # percent of people vaccinated for given strategy
#     if (percent_vax > percent_strat){
#       nvax <- percent_vax*npop - percent_strat*npop
#       leftover_groups <- 1:9
#       leftover_groups <- leftover_groups[!leftover_groups %in% groups]
#       people_to_vax <- sum(N[leftover_groups])
#       vax_proportion <- rep(0, nage)
#       vax_proportion[leftover_groups] <- N[leftover_groups]/people_to_vax
# 
#       vax_distribution <- nvax*vax_proportion
#       vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
# 
#       prob_vaccinated <- vax_distribution/N
# 
#       #V_0 <- V_0 + vax_distribution - N*prob_vaccinated*sero #leaky
#       V_0 <- V_0 + vax_distribution*v_e - N*prob_vaccinated*sero*v_e # all-or-nothing
#     }
#     
#   }
# 
#   # initial I: 1 in each age group unless everyone is vaccinated and/or sero positive
#   I_0 <- rep(0, nage)
#   I_0[(N - R_0 - V_0) > 1] <- 1
#   
#   S_0 <- N - I_0 - V_0 - R_0
#   
#   inits <- c(S_0,E_0,I_0,R_0,V_0)
# 
#   # _____________________________________________________________________
#   # NUMERICALLY SOLVE 
#   # _____________________________________________________________________
#   parameters = list(u=u, d_E=d_E, d_I=d_I, C=C, v_e=v_e)
# 
#   # t <- seq(0,800,1)
#   # tfinal <- 800
#   # 
#   # df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
# 
#   running = TRUE
#   t <- c(0:365)
#   df <- as.data.frame(deSolve::lsoda(inits, t, calculate_derivatives, parameters))
#   ## run til end of outbreak (I < 1 person)
#   # t <- t + 50
#   # 
#   # while(running == TRUE){
#   #   inits <- as.numeric(df[t[1]+1, -(1)])
#   #   temp <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
#   #   row.names(temp) <- t+1
#   #   temp <- temp[-(1),]
#   # 
#   #   df <- rbind(df, temp)
#   # 
#   #   I_tfinal <- sum(df[(t+1),20:28])
#   #   if (I_tfinal < 1){
#   #     running = FALSE
#   #   } else {
#   #     t <- t + 50
#   #     }
#   # }
# 
#   names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
#                     "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
#                     "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
#                     "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
#                     "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
# 
#   return(df)
# }

# run_sim_nontransmissionblocking = function(C, percent_vax, strategy, alpha, omega, u = u_constant, v_e = v_e_constant, frac_age = age_demo, npop = pop_total, N = N_i, nage = num_groups, sero = sero_none, sero_testing = FALSE){
#   
#   # Disease Tranmission
#   d_E    <- 1/3 # incubation period (E -> I), ref: Davies
#   d_I <- 1/5 # recovery period (I -> R), ref: Davies
#   
#   # _____________________________________________________________________
#   # INITIAL CONDITIONS 
#   # Vaccine strategies distribute proportionally to 
#   #     all: all age groups
#   #     kids: age groups < 20
#   #     adults: age groups 20-49
#   #     elderly: age groups 60+
#   #     20+: age groups 20+
#   # _____________________________________________________________________
#   E_0    <- rep(0,nage)
#   R_0    <- N * sero
#   
#   # specify group to vaccinate according to allocation strategy
#   if (strategy == "no vax"){
#     V_0 <- rep(0, nage)
#   } else {
#     if (strategy == "all"){ 
#       groups <- 1:9
#     } else if (strategy == "kids"){ 
#       groups <- 1:2
#     } else if (strategy == "adults") { 
#       groups <- 3:5
#     } else if (strategy == "elderly") {
#       groups <- 7:9
#     } else if (strategy == "20+") {
#       groups <- 3:9
#     }
#     people_to_vax <- sum(N[groups])
#     
#     vax_proportion <- rep(0, nage)
#     vax_proportion[groups] <- N[groups]/people_to_vax
#     
#     nvax <- percent_vax*npop 
#     vax_distribution <- nvax*vax_proportion
#     vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
#     
#     prob_vaccinated <- vax_distribution/N
#     
#     if (sero_testing == FALSE){
#       # account for people who were vaccinated and seropositive
#       #V_0    <- vax_distribution - N*prob_vaccinated*sero
#       V_0 <- vax_distribution*v_e # all or nothing vaccine
#     } else if (sero_testing == TRUE){
#       # assumes everyone vaccinated was seronegative
#       V_0    <- vax_distribution
#       temp <- N - R_0
#       V_0[V_0 > temp] <- temp[V_0 > temp]
#     }
#   }
#   # initial I: 1 in each age group unless everyone is vaccinated and/or sero positive
#   I_0 <- rep(0, nage)
#   I_0[(N - R_0 - V_0) > 1] <- 1
#   
#   S_0 <- N - I_0 - V_0 - R_0
#   
#   inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, SV=V_0, EV = rep(0,nage), IV = rep(0,nage), RV = rep(0,nage))
#   
#   # _____________________________________________________________________
#   # NUMERICALLY SOLVE 
#   # _____________________________________________________________________
#   parameters = list(u=u, d_E = d_E, d_I = d_I, C=C, v_e=v_e, alpha = alpha, omega = omega)
# 
#   t <- seq(0,1000,1)
#   tfinal <- 1000
#   
#   df <- as.data.frame(deSolve::lsoda(inits, t, calculate_derivatives_nontransmissionblocking, parameters))
#   
#   I_tfinal <- sum(df[(tfinal+1),20:28])
#   
#   if (I_tfinal > 1){
#     print("Warning: simulation may not have run long enough")
#   }
#   
#   df
# }

# run_sim_WHOstrat = function(C, percent_vax, V_0, u = u_constant, v_e = v_e_constant, 
#                    frac_age = age_demo, npop = pop_total, N = N_i, nage = num_groups, 
#                    sero = sero_none, sero_testing = FALSE){
#   
#   # Disease Tranmission
#   d_E    <- 1/3 # incubation period (E -> I), ref: Davies
#   d_I <- 1/5 # recovery period (I -> R), ref: Davies
#   
#   E_0    <- rep(0,nage)
#   R_0    <- N * sero
#   
#   # initial I: 1 in each age group unless everyone is vaccinated and/or sero positive
#   I_0 <- rep(0, nage)
#   I_0[(N - R_0 - V_0) > 1] <- 1
#   
#   S_0 <- N - I_0 - V_0 - R_0
#   
#   inits <- c(S_0,E_0,I_0,R_0,V_0)
#   
#   # _____________________________________________________________________
#   # NUMERICALLY SOLVE 
#   # _____________________________________________________________________
#   parameters = list(u=u, d_E=d_E, d_I=d_I, C=C, v_e=v_e)
#   
#   # t <- seq(0,800,1)
#   # tfinal <- 800
#   # 
#   # df <- as.data.frame(lsoda(inits, t, calculate_derivatives, parameters))
#   
#   running = TRUE
#   t <- c(0:50)
#   df <- as.data.frame(deSolve::lsoda(inits, t, calculate_derivatives, parameters))
#   t <- t + 50
#   
#   while(running == TRUE){
#     inits <- as.numeric(df[t[1]+1, -(1)])
#     temp <- as.data.frame(deSolve::lsoda(inits, t, calculate_derivatives, parameters))
#     row.names(temp) <- t+1
#     temp <- temp[-(1),]
#     
#     df <- rbind(df, temp)
#     
#     I_tfinal <- sum(df[(t+1),20:28])
#     if (I_tfinal < 1){
#       running = FALSE
#     } else {
#       t <- t + 50
#     }
#   }
#   
#   names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
#                  "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
#                  "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
#                  "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
#                  "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
#   
#   return(df)
# }
# 
# 
# run_sim_later_vax <- function(t_infected, strategy, percent_vax, sero_testing = FALSE, 
#                               npop = pop_total, N = N_i, nage = num_groups) {
#   # t_infected = time to distribute vaccine wrt % of infections that have occured
#   
#   df_baseline <- list_all[[1]]
#   R_time <- apply(df_baseline[29:37], 1, sum)
#   R_time <- (R_time/pop_total)*100
#   tot_infected <- compute_total_cases(list_all[[1]])
#   t_init <- length(R_time[R_time <= tot_infected*t_infected])
#   
#   inits <- as.numeric(df_baseline[t_init, -(1)])
#   S_0 <- inits[1:9]
#   E_0 <- inits[10:18]
#   I_0 <- inits[19:27]
#   R_0 <- inits[28:36]
#   
#   sero <- R_0/N_i
#   # age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
#   # df <- data.frame(age_groups, sero)
#   # 
#   # theme_set(theme_minimal(base_size = 15))
#   # 
#   #  ggplot(df, aes(y=sero*100, x=age_groups)) +
#   #    geom_bar(position="stack", stat="identity") +
#   #    xlab("Age group") +
#   #    ylab("Percent Seropositive") +
#   #    scale_x_discrete(breaks=c("0-9","20-29", "40-49", "60-69", "80+")) +
#   #    theme(plot.subtitle = element_text(size = 9, face = "italic", hjust = 1),
#   #          plot.title = element_text(hjust = 0.5)) +
#   #    ylim(0, 65)
#   
#           
#   if (strategy == "all"){ 
#     groups <- 1:9
#   } else if (strategy == "kids"){ 
#     groups <- 1:2
#   } else if (strategy == "adults") { 
#     groups <- 3:5
#   } else if (strategy == "elderly") {
#     groups <- 7:9
#   } else if (strategy == "20+") {
#     groups <- 3:9
#   }
#   people_to_vax <- sum(N[groups])
#   
#   vax_proportion <- rep(0, nage)
#   vax_proportion[groups] <- N[groups]/people_to_vax
#   
#   nvax <- percent_vax*npop 
#   vax_distribution <- nvax*vax_proportion
#   vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
#   
#   prob_vaccinated <- vax_distribution/N
#   
#   if (sero_testing == FALSE){
#     # account for people who were vaccinated and seropositive
#     V_0    <- vax_distribution - N*prob_vaccinated*sero
#     V_0[V_0 > S_0] <- S_0[V_0 > S_0]
#   } else if (sero_testing == TRUE){
#     # assumes everyone vaccinated was seronegative
#     V_0    <- vax_distribution
#     V_0[V_0 > S_0] <- S_0[V_0 > S_0]
#   }
#   S_0 <- S_0 - V_0
#   
#   inits <- c(S=S_0,E=E_0,I=I_0,R=R_0, V=V_0)
#   
#   nu    <- 1/3 # incubation period (E -> I), ref: Davies
#   gamma <- 1/5 # recovery period (I -> R), ref: Davies
#   parameters = list(beta=u_constant, nu=nu, gamma=gamma, C=C, v_e=v_e_constant)
# 
#   t <- seq(0,150,1)
# 
#   df <- as.data.frame(deSolve::lsoda(inits, t, calculate_derivatives, parameters))
#   df
# }
# 
# run_sim_vax_overtime = function(C, percent_vax, strategy, num_perday, u = u_constant, v_e = v_e_constant, 
#                    frac_age = age_demo, npop = pop_total, N = N_i, nage = num_groups, 
#                    sero = sero_none, sero_testing = FALSE){
#   # New IC: start w/ 0.5% of each age group in I
#   # Only run sim for 365 days post start of rollout
#   # Vaccine rollout is continuous at 1% of total pop/day until all vaccines are distributed
#   
#   # Disease Tranmission
#   d_E    <- 1/3 # incubation period (E -> I), ref: Davies
#   d_I <- 1/5 # recovery period (I -> R), ref: Davies
#   
#   E_0    <- rep(0,nage)
#   R_0    <- N * sero
#   I_0    <- N*0.005 # 0.5% of each age group
#   
#   # specify group to vaccinate according to allocation strategy
#   if (strategy == "no vax"){
#     V_0 <- rep(0, nage)
#   } else {
#     if (strategy == "all"){ 
#       groups <- 1:9
#     } else if (strategy == "kids"){ 
#       groups <- 1:2
#     } else if (strategy == "adults") { 
#       groups <- 3:5
#     } else if (strategy == "elderly") {
#       groups <- 7:9
#     } else if (strategy == "20+") {
#       groups <- 3:9
#     }
#     people_to_vax <- sum(N[groups])
#     
#     vax_proportion <- rep(0, nage)
#     vax_proportion[groups] <- N[groups]/people_to_vax
#     
#     vax_supply <- percent_vax*npop
# 
#     if (vax_supply > num_perday*npop){
#       nvax <- num_perday*npop
#       vax_supply <- vax_supply - nvax
#     }
#     else {
#       nvax <- vax_supply
#       vax_supply <- 0
#     }
# 
#     vax_distribution <- nvax*vax_proportion
#     vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
# 
#     prob_vaccinated <- vax_distribution/N
# 
#     if (sero_testing == FALSE){
#       # account for people who were vaccinated and seropositive
#       #V_0    <- vax_distribution - N*prob_vaccinated*sero # leaky vaccine
#       V_0 <- vax_distribution*v_e - N*prob_vaccinated*sero # all-or-nothing vaccine
#     } else if (sero_testing == TRUE){
#       # assumes everyone vaccinated was seronegative
#       #V_0    <- vax_distribution # leaky vaccine
#       V_0    <- vax_distribution*v_e # all_or_nothing vaccine
#     }
#       temp <- N - R_0 - I_0
#       V_0[V_0 > temp] <- temp[V_0 > temp]
#   }
#   # nvax_vec = rep(0, 365)
#   # count = 1
#   # while(vax_supply > num_perday*npop){
#   #   nvax_vec[count] <- num_perday*npop
#   #   vax_supply <- vax_supply - num_perday*npop
#   #   count = count + 1
#   # }
#   # nvax_vec[count] <- vax_supply
#   # 
#   # V_0 <- rep(0, nage)
#   # 
#   S_0 <- N - I_0 - V_0 - R_0
#   
#   inits <- c(S_0,E_0,I_0,R_0,V_0)
#   
#   t_switch <- people_to_vax/(npop*num_perday) # number of days it takes to vaccinate everyone in the strategy
#   
#   # _____________________________________________________________________
#   # NUMERICALLY SOLVE 
#   # _____________________________________________________________________
#   parameters = list(u=u, d_E=d_E, d_I=d_I, C=C, v_e=v_e, vax_supply = vax_supply, 
#                     vax_proportion = vax_proportion, npop = npop, N = N, num_perday = num_perday)
#   
#   running = TRUE
#   t <- c(0:1)
#   df <- as.data.frame(deSolve::lsoda(inits, t, calculate_derivatives_vax_overtime, parameters))
#   t <- t + 1
#   
#   vax_supply <- vax_supply - num_perday*npop
#   vax_supply[vax_supply < 0] <- 0
#   
#   while(running == TRUE){
#     inits <- as.numeric(df[t[2], -(1)])
#     
#     if (t[2] == ceiling(t_switch)){
#       leftover_groups <- 1:9
#       people_to_vax <- sum(N[leftover_groups])
#       vax_proportion <- rep(0, nage)
#       vax_proportion[leftover_groups] <- N[leftover_groups]/people_to_vax
#     }
#     else if (t[2] > ceiling(t_switch)){
#       leftover_groups <- 1:9
#       leftover_groups <- leftover_groups[!leftover_groups %in% groups]
#       people_to_vax <- sum(N[leftover_groups])
#       vax_proportion <- rep(0, nage)
#       vax_proportion[leftover_groups] <- N[leftover_groups]/people_to_vax
#     }
#     
#     parameters = list(u=u, d_E=d_E, d_I=d_I, C=C, v_e=v_e, vax_supply=vax_supply, 
#                       vax_proportion = vax_proportion, npop = npop, N = N, num_perday = num_perday)
#     temp <- as.data.frame(deSolve::lsoda(inits, t, calculate_derivatives_vax_overtime, parameters))
#     row.names(temp) <- t+1
#     temp <- temp[-(1),]
#     
#     df <- rbind(df, temp)
#     
#     vax_supply <- vax_supply - num_perday*npop
#     vax_supply[vax_supply < 0] <- 0
#     
#     if (vax_supply == 0){
#       remaining_t = c(t[2]:365)
#       inits <- as.numeric(df[t[2]+1, -(1)])
#       parameters = list(u=u, d_E=d_E, d_I=d_I, C=C, v_e=v_e, vax_supply=vax_supply, 
#                         vax_proportion = vax_proportion, npop = npop, N = N, num_perday = num_perday)
#       temp <- as.data.frame(deSolve::lsoda(inits, remaining_t, calculate_derivatives_vax_overtime, parameters))
#       row.names(temp) <- remaining_t+1
#       temp <- temp[-(1),]
#       
#       df <- rbind(df, temp)
#       running = FALSE
#     } else if (t[2] == 365){
#       running = FALSE
#     } else {
#       t <- t + 1
#     }
#   }
#   
#   names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9",
#                  "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9",
#                  "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9",
#                  "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9",
#                  "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9")
#   
#   return(df)
# }

run_sim_new = function(C, percent_vax, strategy, num_perday, v_e_type, v_e = v_e_constant,
                       u = u_var, sero = sero_none, sp = 1, se = 0, refuse_vax = 0.3, syn_sero_compartments = NA){
  # IC: 0.25% of each age group exposted and infected if numperday != 1
  #     Sero+ start in recovered
  #     refuse_vax % proportion of each compartment start in *x (Sx, Ex, Ix and Rx)
  #     everyone else starts in S
  # Runtime: 365 days post start of rollout
  # Rollout: Continuous at num_perday% of total pop/day until all vaccines are distributed
  #          so num_perday = 1 implies anticipatory rollout
  
  # Disease Tranmission
  d_E <- 1/3 # incubation period (E -> I), ref: Davies
  d_I <- 1/5 # recovery period (I -> R), ref: Davies
  
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
  
  compartments_initial <- c(S_0,Sv_0,Sx_0,E_0,Ev_0,Ex_0,I_0,Iv_0,Ix_0,R_0,Rv_0,Rx_0,D_0)
  
  if (length(syn_sero_compartments) > 1){
    compartments_initial <- syn_sero_compartments
  }
  
  compartments_aftervax <- move_vaccinated(compartments_initial, num_perday, vax_supply, vax_proportion,
                                           groups, v_e, v_e_type, sp, se)
  
  parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_e_type = v_e_type)
  
  running = TRUE
  t <- c(0:1)
  df <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_new, parameters))
  df[1,] <- c(0, compartments_initial)
  
  t <- t + 1
  
  vax_supply <- vax_supply - num_perday*pop_total
  vax_supply[vax_supply < 0] <- 0

  while(running == TRUE){
    compartments_initial <- as.numeric(df[t[2], -(1)])
    
    compartments_aftervax <- move_vaccinated(compartments_initial, num_perday, vax_supply, vax_proportion,
                                             groups, v_e, v_e_type, sp, se)
    
    parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_e_type = v_e_type)
    temp <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_new, parameters))
    row.names(temp) <- t+1
    temp <- temp[-(1),]
    
    df <- rbind(df, temp)
    
    vax_supply <- vax_supply - num_perday*pop_total
    vax_supply[vax_supply < 0] <- 0
    
    if (vax_supply == 0){
      remaining_t = c(t[2]:365)
      inits <- as.numeric(df[t[2]+1, -(1)])
      parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_e_type = v_e_type)
      temp <- as.data.frame(deSolve::lsoda(inits, remaining_t, calculate_derivatives_new, parameters))
      row.names(temp) <- remaining_t+1
      temp <- temp[-(1),]
      
      df <- rbind(df, temp)
      running = FALSE
    } else if (t[2] == 365){
      running = FALSE
    } else {
      t <- t + 1
    }
  }
  
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
                       u = u_var, sero = sero_none, refuse_vax = 0.3){
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
  
  compartments_initial <- c(S_0,Sv_0,Sx_0,E_0,Ev_0,Ex_0,I_0,Iv_0,Ix_0,R_0,Rv_0,Rx_0,D_0)
  
  compartments_aftervax <- move_vaccinated_NTB(compartments_initial, num_perday, vax_supply, 
                                               vax_proportion, groups)
  
  parameters = list(u=u, C=C, ve_S = ve_S, ve_I = ve_I, ve_P = ve_P)
  
  running = TRUE
  t <- c(0:1)
  df <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_NTB, parameters))
  df[1,] <- c(0, compartments_initial)
  
  t <- t + 1
  
  vax_supply <- vax_supply - num_perday*pop_total
  vax_supply[vax_supply < 0] <- 0
  
  while(running == TRUE){
    compartments_initial <- as.numeric(df[t[2], -(1)])
    
    compartments_aftervax <- move_vaccinated_NTB(compartments_initial, num_perday, vax_supply, vax_proportion,
                                             groups)
    
    parameters = list(u=u, C=C, ve_S = ve_S, ve_I = ve_I, ve_P = ve_P)
    temp <- as.data.frame(deSolve::lsoda(compartments_aftervax, t, calculate_derivatives_NTB, parameters))
    row.names(temp) <- t+1
    temp <- temp[-(1),]
    
    df <- rbind(df, temp)
    
    vax_supply <- vax_supply - num_perday*pop_total
    vax_supply[vax_supply < 0] <- 0
    
    if (vax_supply == 0){
      remaining_t = c(t[2]:365)
      inits <- as.numeric(df[t[2]+1, -(1)])
      parameters = list(u=u, C=C, ve_S = ve_S, ve_I = ve_I, ve_P = ve_P)
      temp <- as.data.frame(deSolve::lsoda(inits, remaining_t, calculate_derivatives_NTB, parameters))
      row.names(temp) <- remaining_t+1
      temp <- temp[-(1),]
      
      df <- rbind(df, temp)
      running = FALSE
    } else if (t[2] == 365){
      running = FALSE
    } else {
      t <- t + 1
    }
  }
  
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
