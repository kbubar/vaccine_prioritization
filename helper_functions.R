library(grid)
library(dplyr)
library(data.table)
calculate_derivatives=function(t, x, parameters){
  # the parameters in the parameters list are:
  #    the probability of transmission on contact, beta
  #    the incubation period, nu
  #    the recovery period, gamma
  #    the contact matrix, C, that is the # contacts per day among age groups
  #
  # Note that x is a vector of length (#model compartment types)*(#age classes)
  # Thus, S, E, I and R are vectors, all of length num_groups
  ncompartment <- 5
  num_groups <- length(x)/ncompartment
  S    <- as.matrix(x[1:num_groups])
  E    <- as.matrix(x[(num_groups+1):(2*num_groups)])
  I    <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  R    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  V    <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  
  I[I<0] = 0
 
  u <- parameters$u
  C <- parameters$C
  d_E <- parameters$d_E
  d_I <- parameters$d_I
  v_e <- parameters$v_e
  
  N = S+E+I+R+V
  dS = -as.matrix(S*u)*(as.matrix(C)%*%as.matrix(I/N))
  dV = (1-v_e)*(-as.matrix(V*u))*(as.matrix(C)%*%as.matrix(I/N)) # multiply by 0 if all-or-nothing
  dE = -dS - dV - d_E*as.matrix(E)
  dI = d_E*as.matrix(E) - d_I*as.matrix(I)
  dR = +d_I*as.matrix(I)
  
  out=c(dS,dE,dI,dR,dV)
  list(out)
}

calculate_derivatives_nontransmissionblocking=function(t, x, parameters){
  # the parameters in the parameters list are:
  #    the probability of transmission on contact, beta
  #    the incubation period, nu
  #    the recovery period, gamma
  #    the contact matrix, C, that is the # contacts per day among age groups
  #
  # Note that x is a vector of length (#model compartment types)*(#age classes)
  # Thus, S, E, I and R are vectors, all of length num_groups
  ncompartment <- 8
  num_groups <- length(x)/ncompartment
  S    <- as.matrix(x[1:num_groups])
  E    <- as.matrix(x[(num_groups+1):(2*num_groups)])
  I    <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  R    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  SV    <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  EV    <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  IV    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  RV    <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  
  I[I<0] = 0
  IV[IV<0] = 0
  
  u <- parameters$u
  C <- parameters$C
  d_E <- parameters$d_E
  d_I <- parameters$d_I
  v_e <- parameters$v_e
  alpha <- parameters$alpha
  omega <- parameters$omega
  
  N = S+E+I+R+SV+EV+IV+RV
  dS = -as.matrix(S*u)*(as.matrix(C)%*%as.matrix((I+ (omega*IV))/N))
  dSV = -as.matrix(SV*alpha*u)*(as.matrix(C)%*%as.matrix((I+ (omega*IV))/N))
  dE = -dS - d_E*as.matrix(E)
  dEV = -dSV - d_E*as.matrix(EV)
  dI = d_E*as.matrix(E) - d_I*as.matrix(I)
  dIV = d_E*as.matrix(EV) - d_I*as.matrix(IV)
  dR = +d_I*as.matrix(I)
  dRV = +d_I*as.matrix(IV)
  
  out=c(dS,dE,dI,dR,dSV,dEV,dIV,dRV)
  list(out)
}

get_v_e = function(p, y0, hinge_age){
  # INPUT: p is the final v_e % (as a decimal) 
  #       i.e. p = 1 -> perfect vaccine at all ages
  #       y0 is the baseline v_e 
  #       hinge_age is the age group that v_e begins decreasing after
  # ASSUMPTIONS: perfect v_e up to hinge age then linear decline to p
  # OUTPUT: vector v_e 
  
  x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
  y <- y0 + ((p-y0)/(80-hinge_age))*(x - hinge_age)
  index <- match(c(hinge_age), x)
  y[1:index] <- y0
  
  y
}

get_best_strat = function(all, kids, adults, 
                          elderly, twentyplus){
  
  total_deaths <- c(compute_total_deaths(all$`15`),
                    compute_total_deaths(kids$`15`), 
                    compute_total_deaths(adults$`15`),
                    compute_total_deaths(elderly$`15`), 
                    compute_total_deaths(twentyplus$`15`))

  baseline_deaths <- compute_total_deaths(all$`0`)
  reduction_in_deaths <- (1 - (total_deaths/baseline_deaths))*100
  
  strats <- c("all", "kids", "adults", "elderly", "twentyplus")
  strats[which.max(reduction_in_deaths)]
}

run_v_e_var = function(p, v_e_baseline, hinge_age){
  v_e_var <- get_v_e(p, y0 = v_e_baseline, hinge_age = hinge_age)
  
  all <- vector(mode = "list")
  kids <- vector(mode = "list")
  adults <- vector(mode = "list")
  elderly <- vector(mode = "list")
  twentyplus <- vector(mode = "list")
  
  for (i in seq(0, 15, by = 15)){
    j <- i/100
    all[[paste0(i)]] <- run_sim(C, j, "all", u_var, v_e_var)
    kids[[paste0(i)]] <- run_sim(C, j, "kids", u_var, v_e_var)
    adults[[paste0(i)]] <- run_sim(C, j, "adults", u_var, v_e_var)
    elderly[[paste0(i)]] <- run_sim(C, j, "elderly", u_var, v_e_var)
    twentyplus[[paste0(i)]] <- run_sim(C, j, "20+", u_var, v_e_var)
  }
  
  return(get_best_strat(all, kids, adults, elderly, twentyplus))
}

v_e_bisection = function(v_e_baseline, hinge_age){
  p_high <- 0.5
  p_low <- 0
  
  temp <- run_v_e_var(p_high, v_e_baseline, hinge_age = hinge_age)
  if (temp != "elderly"){
    print("Warning: elderly is not most strategic to start")
    return(0.5)
  }
  
  temp <- run_v_e_var(p_low, v_e_baseline, hinge_age = hinge_age)
  if (temp == "elderly"){
    print("Warning: elderly is always most strategic")
    return(0)
  }
  
  running = TRUE
  while(running){
    p_next = (p_high + p_low)/2
    temp <- run_v_e_var(p_next, v_e_baseline, hinge_age)
    if (temp == "elderly"){
      p_high <- p_next
    } else {
      p_low <- p_next
    }
    
    if ((p_high - p_low) < 0.01){
      return(p_next)
    }
  }
  
}

compute_R0 = function(u, C){
  gamma <- 1/5 # recovery period (I -> R), ref: Davies
  
  # Davies NGM
  Du <- diag(u, 9)
  Dy <- diag(1/gamma, 9)
  
  NGM <- Du %*% C %*% Dy
  R0  <- abs(eigen(NGM)$values[1])
}

compute_total_cases = function(df){
  infections <- rep(0,num_groups) # number of age groups
  
  R_index <- 29 # col number for R1
  
  for (i in 1:num_groups) {
    # infections = total # recovered - initial recovered (seropositive)
    infections[i] <- max(df[R_index]) - min(df[R_index])
    R_index <- R_index + 1
  }
  
  tot_infections <- sum(infections)/pop_total * 100
  #tot_infections <- round(tot_infections, 1)
}

compute_total_deaths = function(df){
  deaths <- rep(0,num_groups)
  
  R_index <- 29 # col number for R1
  
  for (i in 1:num_groups) {
    # infections = total # recovered - initial recovered (seropositive)
    deaths[i] <- (max(df[R_index]) - min(df[R_index]))*IFR[i]
    R_index <- R_index + 1
  }
  
  tot_deaths <- sum(deaths)/pop_total * 100
}

compute_total_YLL = function(df){
  YLL <- rep(0,num_groups)
  
  R_index <- 29 # col number for R1
  
  for (i in 1:num_groups) {
    # infections = total # recovered - initial recovered (seropositive)
    deaths_temp <- (max(df[R_index]) - min(df[R_index]))*IFR[i]
    YLL[i] <- YLL_vec[i]*deaths_temp
    R_index <- R_index + 1
  }
  
  tot_YLL <- sum(YLL) # absolute number of YLL
}


get_reduction_in_cases_df = function(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  total_cases <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  
  num_strategies <- 5
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies*2))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list*2), rep("kids", num_per_list*2), 
             rep("adults", num_per_list*2), rep("elderly", num_per_list*2), 
             rep("twentyplus", num_per_list*2))
  temp <-  c(rep("constant", num_per_list), rep("var", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  baseline_cases <- compute_total_cases(list_all$`0`)
  baseline_cases_var <- compute_total_cases(list_all_var$`0`)
  
  temp <- c(rep(baseline_cases, num_per_list), rep(baseline_cases_var, num_per_list))
  baseline_cases <- c(rep(temp, num_strategies))
  
  reduction <- (1-(total_cases/baseline_cases))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

get_reduction_in_cases_df_novar = function(){
  total_cases <- rep(NA, 255)
  
  count <- 1
  for (i in list_all){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_cases[count] <- compute_total_cases(i)
    count <- count + 1
  }
  
  num_strategies <- 5
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list), rep("kids", num_per_list), 
             rep("adults", num_per_list), rep("elderly", num_per_list), 
             rep("twentyplus", num_per_list))
  temp <-  c(rep("constant", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  baseline_cases <- compute_total_cases(list_all$`0`)
  
  temp <- c(rep(baseline_cases, num_per_list))
  baseline_cases <- c(rep(temp, num_strategies))
  
  reduction <- (1-(total_cases/baseline_cases))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

get_reduction_in_deaths_df = function(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  total_deaths <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  
  num_strategies <- 5
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies*2))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list*2), rep("kids", num_per_list*2), 
             rep("adults", num_per_list*2), rep("elderly", num_per_list*2), 
             rep("twentyplus", num_per_list*2))
  temp <-  c(rep("constant", num_per_list), rep("var", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  baseline_deaths <- compute_total_deaths(list_all$`0`)
  baseline_deaths_var <- compute_total_deaths(list_all_var$`0`)
  
  temp <- c(rep(baseline_deaths, num_per_list), rep(baseline_deaths_var, num_per_list))
  baseline_deaths <- c(rep(temp, num_strategies))
  
  reduction <- (1 - (total_deaths/baseline_deaths))*100
  
  df <- data.frame(vax_avail, strat, total_deaths, reduction, variable)
}

get_reduction_in_deaths_df_novar = function(){
  total_deaths <- rep(NA, 255)
  
  count <- 1
  for (i in list_all){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  
  num_strategies <- 5
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list), rep("kids", num_per_list), 
             rep("adults", num_per_list), rep("elderly", num_per_list), 
             rep("twentyplus", num_per_list))
  temp <-  c(rep("constant", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  baseline_deaths <- compute_total_deaths(list_all$`0`)
  
  temp <- c(rep(baseline_deaths, num_per_list))
  baseline_deaths <- c(rep(temp, num_strategies))
  
  reduction_in_deaths <- (1 - (total_deaths/baseline_deaths))*100
  
  df <- data.frame(vax_avail, strat, total_deaths, reduction_in_deaths, variable)
}

get_reduction_in_YLL_df = function(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  total_YLL <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total_YLL[count] <- compute_total_YLL(i)
    count <- count + 1
  }
  
  num_strategies <- 5
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies*2))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list*2), rep("kids", num_per_list*2), 
             rep("adults", num_per_list*2), rep("elderly", num_per_list*2), 
             rep("twentyplus", num_per_list*2))
  temp <-  c(rep("constant", num_per_list), rep("var", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  baseline_YLL <- compute_total_YLL(list_all$`0`)
  baseline_YLL_var <- compute_total_YLL(list_all_var$`0`)
  
  temp <- c(rep(baseline_YLL, num_per_list), rep(baseline_YLL_var, num_per_list))
  baseline_YLL <- c(rep(temp, num_strategies))
  
  reduction <- (1 - (total_YLL/baseline_YLL))*100
  
  df <- data.frame(vax_avail, strat, total_YLL, reduction, variable)
}
# Plot trajectories ----
plot_allages_onestrategy =function(df, compartment){
  # INPUTS:
  # df: simulation to plot
  # compartment: character of compartment to plot i.e. "I", "R"
  #
  # OUTPUT:
  # p: ggplot object of plot
  
  # Dataframe for total number of infections/recoveries
  final_time <- as.numeric(dim(df)[1])-1
  
  total_df <- data.frame(time = 0:final_time,
                      Infected = apply(df[20:28], 1, sum),
                      Recovered = apply(df[29:37], 1, sum))

  total_df <- total_df %>%
    select(time, Infected) %>%
    gather(key = "Compartment", value = "num", -time)
  
  total_df$percent <- total_df$num/pop_total*100
  
  # Dataframe for infections by age
  col_to_gather <- vector(mode="character", length=num_groups)
  for (i in 1:9) {col_to_gather[i] <- paste0(compartment,i)}

  new_df <- df %>%
    select(time, all_of(col_to_gather)) %>%
    gather(key = "age_group", value = "num", -time)

  # Add col for percent infected & recovered
  new_df$percent <- new_df$num
  for (i in 1:length(new_df$num)){
    if (new_df$age_group[i] == col_to_gather[1]) {new_df$percent[i] <- new_df$percent[i]/N[1]}
    else if (new_df$age_group[i] == col_to_gather[2]) {new_df$percent[i] <- new_df$percent[i]/N[2]}
    else if (new_df$age_group[i] == col_to_gather[3]) {new_df$percent[i] <- new_df$percent[i]/N[3]}
    else if (new_df$age_group[i] == col_to_gather[4]) {new_df$percent[i] <- new_df$percent[i]/N[4]}
    else if (new_df$age_group[i] == col_to_gather[5]) {new_df$percent[i] <- new_df$percent[i]/N[5]}
    else if (new_df$age_group[i] == col_to_gather[6]) {new_df$percent[i] <- new_df$percent[i]/N[6]}
    else if (new_df$age_group[i] == col_to_gather[7]) {new_df$percent[i] <- new_df$percent[i]/N[7]}
    else if (new_df$age_group[i] == col_to_gather[8]) {new_df$percent[i] <- new_df$percent[i]/N[8]}
    else if (new_df$age_group[i] == col_to_gather[9]) {new_df$percent[i] <- new_df$percent[i]/N[9]}
  }
  new_df$percent <- new_df$percent * 100

  # Plot
  theme_set(theme_minimal(base_size = 20))
  
  string_to_add <- paste0("Total infected: ", round(compute_total_cases(df),3), "%")
  
  # Create a text
  grob <- grobTree(textGrob(string_to_add, x=0.5,  y=0.95, hjust=0,
                            gp=gpar(col="white", fontsize=13, fontface="bold")))

  # p <- ggplot(total_df, aes(x = time, y = percent)) +
  #   theme_classic(base_size = 25) +
  #   geom_hline(yintercept = max(total_df$percent), color = "#b7b7b7ff", linetype = "dashed", size = 1.7) +
  #   geom_line(color = "#595959", size = 2) +
  #   xlab("Time") +
  #   scale_x_continuous(expand = c(0,0),limit = c(0,400), breaks = c(0,40,80,120)) +
  #   scale_y_continuous(expand = c(0,0),limit = c(0,20), breaks=c(0, round(max(total_df$percent),1), 20)) +
  #   ylab("Infected (%)")

  p <- ggplot(new_df, aes(x = time, y = percent)) +
    theme_classic(base_size = 25) +
    geom_line(aes(color = age_group), size = 2) +
    xlab("Time") +
    scale_x_continuous(expand = c(0,0),limit = c(0,final_time)) +
    # annotation_custom(grob) +
    scale_color_brewer(palette = "Spectral", name = "Age Group",
                       labels =  c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))

  if (compartment == "I") {
    p <- p + ylab("Infected (%)") +
      scale_y_continuous(expand = c(0,0), limit = c(0,25), breaks=c(0, 10, 20))
  } else if (compartment == "R") {
    p <- p + ylab("Percent Recovered") + ylim(0,100)
  }
}

plot_totalinfections_sero =function(df_baseline, df_sero_notest, df_sero_test, compartment){
  # INPUTS:
  # df: simulation to plot
  # compartment: character of compartment to plot i.e. "I", "R"
  #
  # OUTPUT:
  # p: ggplot object of plot
  
  # make new df for plotting
  if (compartment == "I"){
    total <- data.frame(time = 0:150,
                        baseline = apply(df_baseline[20:28], 1, sum), 
                        sero_notest = apply(df_sero_notest[20:28], 1, sum),
                        sero_test = apply(df_sero_test[20:28], 1, sum))
  } else if (compartment == "R"){
    total <- data.frame(time = 0:150,
                        baseline = apply(df_baseline[29:37], 1, sum), 
                        sero_notest = apply(df_sero_notest[29:37], 1, sum),
                        sero_test = apply(df_sero_test[29:37], 1, sum))
  }
  
  new_df <- total %>%
    select(time, baseline, sero_notest, sero_test) %>%
    gather(key = "scenario", value = "num", -time)
  
  # Add col for percent infected & recovered
  new_df$percent <- new_df$num/pop_total*100
  
  # Plot
  theme_set(theme_minimal(base_size = 20))
  
  # string_to_add <- paste0("Total infected: ", compute_total_cases(df), "%")
  
  # Create a text
  # grob <- grobTree(textGrob(string_to_add, x=0.5,  y=0.95, hjust=0,
  #                           gp=gpar(col="white", fontsize=13, fontface="bold")))
  
  # Plot
  p <- ggplot(new_df, aes(x = time, y = percent)) + 
    geom_line(aes(color = scenario), size = 2.5) +
    xlab("Time") +
    guides(colour = guide_legend(override.aes = list(size=6))) +
  #   annotation_custom(grob) +
    scale_color_manual(values = wes_palette("GrandBudapest1", n = 3), name = "Scenario", 
                        labels =  c("Baseline\n(no seropositive)\n", "Seropositive\nno testing\n", "Seropositive\nwith testing"))
  # 
  if (compartment == "I") {
    p <- p + ylab("Percent Infected") + ylim(0, 16)
  } else if (compartment == "R") {
    p <- p + ylab("Percent Recovered") + ylim(0,1)
  }
}

plot_oneage_allstrategy = function(col_name, age_group_num) {
  # INPUTS:
  # col_name: "I2" == infected, age group 2
  # age_group_num: age group of interest, 1 = 0-9, 2 = 10-18, ...
  #
  # OUTPUT:
  # p: ggplot object plotting all strats for one age group & one compartment
  
  nstrat <- 3
  new_df <- data.frame(matrix(ncol = nstrat+1, nrow = dim(df_novax)[1]))
  colnames(new_df) <- c("time", "no_vax", "everyone", "kids")
  new_df$time <- df_novax$time
  
  # Make sure col of interest are stored as vectors
  new_df$no_vax <- unlist(as.data.frame(df_novax)[col_name]/N[age_group_num])
  new_df$everyone <- unlist(as.data.frame(df_prop)[col_name]/N[age_group_num])
  new_df$kids <- unlist(as.data.frame(df_propkids)[col_name]/N[age_group_num])
  
  new_df <- new_df %>%
    select(time, no_vax, everyone, kids) %>%
    gather(key = "strat", value = "percent", -time)
  
  theme_set(theme_minimal(base_size = 15))
  
  p <- ggplot(new_df, aes(x = time, y = percent)) + 
    geom_line(aes(color = strat), size = 2) +
    xlab("Time") +
    scale_color_brewer(palette = "Dark2", name = "Scenario", 
                       labels =  c("Everyone", "Just Kids", "No vaccines"))
  
  if (substr(col_name, 1, 1) == "I"){
    p <- p + ylab("Percent Infected") + ylim(0, 0.4)
  }
}

# Bar plots ----
barplot_vax_strat = function(strategy){
  if (strategy == "no vax"){
    vaccinated <- rep(0,num_groups) * 100
    this_color <- "#808080"
    plot_title <- "No Vaccines"
  } else {
    if (strategy == "all"){
      groups <- 1:9
      this_color <- col_all
      plot_title <- "All Ages"
    } 
    else if (strategy == "kids"){
      groups <- 1:2
      this_color <- col_kids
      plot_title <- "Under 20"
    } 
    else if (strategy == "adults"){
      groups <- 3:5
      this_color <- col_youngadults
      plot_title <- "Adults 20-49"
    } 
    else if (strategy == "elderly"){
      groups <- 7:9
      this_color <- col_elderly
      plot_title <- "Adults 60+"
    }
    else if (strategy == "20+"){
      groups <- 3:9
      this_color <- col_adults
      plot_title <- "Adults 20+"
    }
    people_to_vax <- sum(N_i[groups])
    
    vax_proportion <- rep(0, num_groups)
    vax_proportion[groups] <- N_i[groups]/people_to_vax
  }
  age_groups <- c("0","10", "20", "30", "40", "50", "60", "70", "80")
  vax_proportion <- vax_proportion*100
  
  df <- data.frame(age_groups, vax_proportion, age_demo)
  
  theme_set(theme_classic(base_size = 22))
  
  # plot
  p <- ggplot(df, aes(x=age_groups)) + 
    #geom_bar(aes(y=age_demo), fill="white",color="black", position="stack", stat="identity") + 
    geom_bar(aes(y=vax_proportion), position="stack", stat="identity", fill = this_color, alpha = 0.8) + 
    xlab("Age group") +
    ggtitle(plot_title) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 51), breaks = c(0, 25, 50))
    
  if (strategy == "all"){
    p <- p + scale_x_discrete(breaks=c("0","10", "20", "30", "40", "50", "60", "70", "80"), 
                              labels=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")) +
              theme(plot.title = element_text(color = this_color, size = 19, face = "bold"),
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(angle = 45, hjust = 1),
                      axis.title.y = element_blank())
  } else {
    p <- p + theme(plot.title = element_text(color = this_color,  size = 19, face = "bold"),
                   axis.title.x =element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.y = element_blank()) + 
      scale_x_discrete(breaks=c("0","10", "20", "30", "40", "50", "60", "70", "80"), 
                       labels=c("","","","","","","","",""))  
  }
  return(p)
}

barplot_at_finalT = function(solution_df, IFR, strategy, y, nvax = 0){
  if (strategy == "no vax"){
    vaccinated <- rep(0,num_groups) * 100
    this_color <- "#808080"
    plot_title <- "No Vaccines"
  } 
  else {
    if (strategy == "all"){
      groups <- 1:9
      this_color <- col_all
      plot_title <- "All Ages"
      sub_title <- ""
    } 
    else if (strategy == "kids"){
      groups <- 1:2
      this_color <- col_kids
      plot_title <- "Children (0-19)"
      sub_title <- ""
    } 
    else if (strategy == "adults"){
      groups <- 3:5
      this_color <- col_youngadults
      plot_title <- "Young Adults (20-49)"
      sub_title <- "(20-49)"
    } 
    else if (strategy == "elderly"){
      groups <- 7:9
      this_color <- col_elderly
      plot_title <- "Elderly (60+)"
      sub_title <- "(60+)"
    }
    else if (strategy == "20+"){
      groups <- 3:9
      this_color <- col_adults
      plot_title <- "Adults (20+)"
      sub_title <- "(20+)"
    }
    people_to_vax <- sum(N_i[groups])
    
    vax_proportion <- rep(0, num_groups)
    vax_proportion[groups] <- N_i[groups]/people_to_vax
    
    vax_distribution <- nvax*vax_proportion
    vax_distribution[vax_distribution > N_i] <- N_i[vax_distribution > N_i]
    vaccinated <- vax_distribution/N_i * 100
  }
  
  infections <- rep(0,num_groups)
  deaths <- rep(0, num_groups)
  final_S <- rep(0, num_groups)
  
  R_index <- 29 # col number for R1

  for (i in 1:num_groups) {
    deaths[i] <- max(solution_df[R_index])*IFR[i]/N_i[i]*100
    infections[i] <- max(solution_df[R_index]) # total infected (including dead)
    # infections[i] <- max(solution_df[R_index])*(1-IFR[i])/N_i[i] # total infected (excluding dead)
    final_S[i] <- solution_df[(max(t) + 1), (i + 1)]/N_i[i]
    R_index <- R_index + 1
  }
  
  #age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  age_groups <- c("0","10", "20", "30", "40", "50", "60", "70", "80")
  
  df <- data.frame(age_groups, final_S, deaths, infections, vaccinated)
  
  new_df <- df %>%
    # select(age_groups, infections, deaths, vaccinated, final_S) %>%
    select(age_groups, all_of(y)) %>%
    gather(key = "type", value = "num", -age_groups)
  
  theme_set(theme_classic(base_size = 20))
  
  # plot
  p <- ggplot(new_df, aes(fill=type, y=num, x=age_groups)) + 
    geom_bar(position="stack", stat="identity", fill = this_color) + 
    xlab("Age group") +
    scale_x_discrete(breaks=c("0","10", "20", "30", "40", "50", "60","70", "80"), 
                     labels=c("0","", "20", "", "40", "", "60","", "80")) + 
    scale_y_continuous(expand = c(0,0), limit = c(0, 100), breaks = c(0, 50, 100)) +
    theme(plot.subtitle = element_text(size = 9, face = "italic", hjust = 1),
          plot.title = element_text(size = 17),
          axis.title.x =element_blank(), 
          axis.title.y = element_blank())
  
  if (y == "infections") {
      string_to_add <- paste0("Total infected: ", round(compute_total_cases(solution_df),2), "%")

      p <- p + 
        #ylim(0, 100) + 
        ylab("Infected (%)") + 
        # labs(title = plot_title, subtitle = string_to_add)
        labs(subtitle = string_to_add)
    } 
  else if (y == "deaths") {
      string_to_add <- paste0("Total deaths: ", round(compute_total_deaths(solution_df),2), "%")
      
      p <- p + 
        ylim(0, 8) + 
        ylab("Deaths (%)") + 
        # labs(title = plot_title, subtitle = string_to_add)
        labs(subtitle = string_to_add)
    } 
  else if (y == "vaccinated"){
      p <- p + 
        #ylim(0, 100) + 
        ylab("Vaccinated (%)") 
        #labs(title = plot_title)
  }
  p # output ggplot
}

barplot_totalcases = function() {
  novax <- compute_total_cases(df_novax)
  everyone_05 <- compute_total_cases(df_prop_05)
  kids_05 <- compute_total_cases(df_propkids_05)
  adults_05 <- compute_total_cases(df_propadults_05)
  everyone_10 <- compute_total_cases(df_prop_10)
  kids_10 <- compute_total_cases(df_propkids_10)
  adults_10 <- compute_total_cases(df_propadults_10)
  
  tot_cases <- c(novax, everyone_05, kids_05, adults_05, everyone_10, kids_10, everyone_10)
  groups <- c("No vax", "Everyone", "Kids", "Adults", "Everyone", "Kids", "Adults")
  vax_avaliable <- c("None", "5%", "5%", "5%", "10%", "10%", "10%")
  df <- data.frame(groups, tot_cases, vax_avaliable)
  
  ggplot(df, aes(fill = vax_avaliable, x=groups, y=tot_cases)) + 
       geom_bar(position = "dodge", stat = "identity") + 
    xlab("Strategy") + 
    ylab("Total infected (%)")
}

barplot_totaldeaths = function() {
  novax <- compute_total_deaths(df_novax)
  everyone_05 <- compute_total_deaths(df_prop_05)
  kids_05 <- compute_total_deaths(df_propkids_05)
  adults_05 <- compute_total_deaths(df_propadults_05)
  everyone_10 <- compute_total_deaths(df_prop_10)
  kids_10 <- compute_total_deaths(df_propkids_10)
  adults_10 <- compute_total_deaths(df_propadults_10)
  
  tot_deaths <- c(novax, everyone_05, kids_05, adults_05, everyone_10, kids_10, everyone_10)
  groups <- c("No vax", "Everyone", "Kids", "Adults", "Everyone", "Kids", "Adults")
  vax_avaliable <- c("None", "5%", "5%", "5%", "10%", "10%", "10%")
  df <- data.frame(groups, tot_deaths, vax)
  
  ggplot(df, aes(fill = vax_avaliable, x=groups, y=tot_deaths)) + 
    geom_bar(position = "dodge", stat = "identity") + 
    xlab("Strategy") + 
    ylab("Total deaths (%)")
}

# Plot over vaccine avaliablity ----
plot_over_vax_avail = function(outcome, var_name, list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
 
  theme_set(theme_minimal(base_size = 26))
  #theme_set(theme_bw(base_size = 26))
  
  percent_kids <- round(sum(age_demo[1:2])*100)
  percent_elderly <- round(sum(age_demo[7:9])*100)
  percent_adults <- round(sum(age_demo[3:5])*100)
  percent_twentyplus <- round(sum(age_demo[3:9])*100)
  
  if (outcome == "cases"){
    df <- get_reduction_in_cases_df(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)
    } 
  else if (outcome == "deaths"){
    df <- get_reduction_in_deaths_df(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)
    #df <- rbind(df[df$variable == "var", ], optimal_df_C[optimal_df_C$variable == "var", ])
  }
  else if (outcome == "YLL"){
    df <- get_reduction_in_YLL_df(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)
  }

  if (percent_kids < 50){df[df$strat == "kids" & df$vax_avail > percent_kids, ]$reduction = NA}
  if (percent_elderly < 50) {df[df$strat == "elderly" & df$vax_avail > percent_elderly, ]$reduction = NA}
  if (percent_adults < 50) {df[df$strat == "adults" & df$vax_avail > percent_adults, ]$reduction = NA}
  if (percent_twentyplus < 50) {df[df$strat == "twentyplus" & df$vax_avail > percent_twentyplus, ]$reduction = NA}
  
  p <- ggplot(df, aes(x = vax_avail, y = reduction, col = strat, fill = strat)) +
    geom_line(size = 2, alpha = 1) + #aes(linetype = variable)
    xlab("Total vaccine supply (% of pop)") +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Adults 20-49", "All Ages", "Adults 60+", 
                                   "Under 20", "Adults 20+")) +
    scale_linetype_discrete(name = var_name,
                            labels = c("Constant", "Age-dependent")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100.2)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    # geom_abline(slope = 1, intercept = 0, size = 3, alpha = 0.2) +
    ggtitle(paste0(country)) +
    theme(legend.position = "none",
          #       axis.title.x =element_blank(),
          #       axis.text.x  = element_blank(),
          #       axis.text.y  = element_blank(),
          #       axis.title.y = element_blank()
    ) +
    theme(legend.text = element_text(size=22),
          legend.title = element_text(size =22),
          title = element_text(size = 22)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
  
  if (outcome == "cases"){ p <- p + ylab("Reduction in infections (%)")}
  else if (outcome == "deaths"){ p <- p + ylab("Reduction in deaths (%)")}
  else if (outcome == "YLL"){ p <- p + ylab("Reduction in YLL (%)")}
  plot(p)
}

plot_over_vax_avail_varyingR0 = function(list_df, outcome, this_R0, this_panel){
  R0_temp <- 2.1
  for (i in 1:3){
    R0_char <- as.character(R0_temp)
    list_df[[i]]$R0 <- R0_char
    R0_temp <- R0_temp + 0.5
  }
  
  df <- do.call(rbind, list_df)
  if (outcome == "cases") {colnames(df)[3] <- "reduction"
  } else {colnames(df)[4] <- "reduction"}
  
  percent_kids <- round(sum(age_demo[1:2])*100)
  percent_elderly <- round(sum(age_demo[7:9])*100)
  percent_adults <- round(sum(age_demo[3:5])*100)
  percent_twentyplus <- round(sum(age_demo[3:9])*100)
  
  if (percent_kids < 50){df[df$strat == "kids" & df$vax_avail > percent_kids, ]$reduction = NA}
  if (percent_elderly < 50) {df[df$strat == "elderly" & df$vax_avail > percent_elderly, ]$reduction = NA}
  if (percent_adults < 50) {df[df$strat == "adults" & df$vax_avail > percent_adults, ]$reduction = NA}
  if (percent_twentyplus < 50) {df[df$strat == "twentyplus" & df$vax_avail > percent_twentyplus, ]$reduction = NA}
  
  theme_set(theme_minimal(base_size = 12))
  
  p <- ggplot() +
    xlab("Total vaccine supply (% of pop)") +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Young Adults (20-49)", "All Ages", "Elderly (60+)", 
                                   "Children (0-19)", "Adults (20+)", "Optimal")) +
    scale_fill_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Young Adults (20-49)", "All Ages", "Elderly (60+)", 
                                   "Children (0-19)", "Adults (20+)", "Optimal")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 101)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
  
    if (this_panel == 1){
      p <- p + theme(axis.title.x = element_blank(),         # 1,2,3,4,5,6
                    #axis.title.y = element_text(size = 20),  # 1,4
                    axis.title.y = element_blank(),
                    axis.text.x = element_blank(),           # 1,2,3
                    axis.text.y = element_text(size = 17),   # 1,4
                    legend.position = "none") 
        #ylab("Reduction in deaths (%)")
    } else if (this_panel == 2 | this_panel == 3){
      p <- p + theme(axis.title.x = element_blank(),         # 1,2,3,4,5,6
                     axis.title.y = element_blank(),         # 2,3,5,6
                    axis.text.x = element_blank(),           # 1,2,3
                    axis.text.y = element_blank(),           # 2,3,5,6
                    legend.position = "none")  
        #ylab("Reduction in deaths (%)")
    } else if (this_panel == 4) {
      p <- p + theme(axis.title.x = element_blank(),          # 1,2,3,4,5,6
                     #axis.title.y = element_text(size = 20),  # 1,4
                     axis.title.y = element_blank(),
                     axis.text.x = element_text(size = 17),   # 4,5,6
                     axis.text.y = element_text(size = 17),   # 1,4
                     legend.position = "none")  
        #ylab("Reduction in infections (%)")
    } else if (this_panel == 5 | this_panel == 6){
      p <- p + theme(axis.title.x = element_blank(),          # 1,2,3,4,5,6
                     axis.title.y = element_blank(),          # 2,3,5,6
                     axis.text.x = element_text(size = 17),   # 4,5,6
                     axis.text.y = element_blank(),           # 2,3,5,6
                     legend.position = "none")   
        #ylab("Reduction in infections (%)")

    }

  strategies <- c("kids","adults", "elderly", "twentyplus", "all")
  j <- 1
  # for (ii in strategies){
  #   p <- p + geom_ribbon(data = df[(df$R0 == "2.1"),],
  #                      aes(x = vax_avail,
  #                          ymin = df[df$R0 == "2.1",]$reduction,
  #                          ymax = df[df$R0 == "3.1",]$reduction,
  #                          fill = strat), 
  #                      alpha = 0.05 
  #                      )
  #   j <- j + 1
  # }
  # for (i in seq(2.1,3.1, by = 0.5)){
  #   p <- p + geom_line(data = df[df$R0 == paste0(i),], 
  #                      aes(x = vax_avail, y = reduction, col = strat), 
  #                      size = 1.1, alpha = 0.15)
  # }
  p <- p + geom_line(data = df[df$R0 == this_R0,], 
                     aes(x = vax_avail, y = reduction, col = strat),
                     size = 1.4, alpha = 1)
  return(p)
} 

plot_over_vax_avail_varyingprev = function(list_df, outcome){
  prev_temp <- 1
  for (i in 1:100){
    prev_char <- as.character(prev_temp)
    list_df[[i]]$prev_num <- prev_char
    prev_temp <- prev_temp + 0.1
  }
  
  df <- do.call(rbind, list_df)
  if (outcome == "cases") {colnames(df)[3] <- "reduction"
  } else {colnames(df)[4] <- "reduction"}
  
  percent_kids <- round(sum(age_demo[1:2])*100)
  percent_elderly <- round(sum(age_demo[7:9])*100)
  percent_adults <- round(sum(age_demo[3:5])*100)
  percent_twentyplus <- round(sum(age_demo[3:9])*100)
  
  if (percent_kids < 50){df[df$strat == "kids" & df$vax_avail > percent_kids, ]$reduction = NA}
  if (percent_elderly < 50) {df[df$strat == "elderly" & df$vax_avail > percent_elderly, ]$reduction = NA}
  if (percent_adults < 50) {df[df$strat == "adults" & df$vax_avail > percent_adults, ]$reduction = NA}
  if (percent_twentyplus < 50) {df[df$strat == "twentyplus" & df$vax_avail > percent_twentyplus, ]$reduction = NA}
  
  theme_set(theme_minimal(base_size = 12))
  
  p <- ggplot() +
    xlab("Total vaccine supply (% of pop)") +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Young Adults (20-49)", "All Ages", "Elderly (60+)", 
                                   "Children (0-19)", "Adults (20+)", "Optimal")) +
    scale_fill_brewer(palette = "Dark2", name = "Allocation Strategy",
                      labels =  c("Young Adults (20-49)", "All Ages", "Elderly (60+)", 
                                  "Children (0-19)", "Adults (20+)", "Optimal")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100), breaks = c(0,25,50, 75, 100)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50), breaks = c(0,10, 20, 30, 40, 50)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
    # ggtitle("")
    
  if (outcome == "cases"){
    p <- p + theme(axis.title.x = element_blank(),         
                   axis.title.y = element_blank(),
                   axis.text.x = element_text(size = 17), 
                   axis.text.y = element_text(size = 17),
                   plot.title = element_blank(),
                   legend.position = "none") 
  } else {
    p <- p + theme(axis.title.x = element_blank(),         
                   axis.title.y = element_blank(), 
                   axis.text.x = element_blank(), 
                   axis.text.y = element_text(size = 17),
                   plot.title = element_text(size = 17),
                   legend.position = "none") 
  }

  for (i in 1:100){
    p <- p + geom_line(data = df[df$prev_num == paste0(i),],
                       aes(x = vax_avail, y = reduction, col = strat),
                       size = 1.1, alpha = 0.25)
  }
  
  if (outcome == "cases") {
    df_Herzog <- get_reduction_in_cases_df(list_all_sero_notest, list_kids_sero_notest, list_adults_sero_notest, list_elderly_sero_notest, list_twentyplus_sero_notest)
    colnames(df_Herzog)[3] <- "reduction"
  } else {
    df_Herzog <- get_reduction_in_deaths_df(list_all_sero_notest, list_kids_sero_notest, list_adults_sero_notest, list_elderly_sero_notest, list_twentyplus_sero_notest)
    colnames(df_Herzog)[4] <- "reduction"}
  
  if (percent_kids < 50){df_Herzog[df_Herzog$strat == "kids" & df_Herzog$vax_avail > percent_kids, ]$reduction = NA}
  if (percent_elderly < 50) {df_Herzog[df_Herzog$strat == "elderly" & df_Herzog$vax_avail > percent_elderly, ]$reduction = NA}
  if (percent_adults < 50) {df_Herzog[df_Herzog$strat == "adults" & df_Herzog$vax_avail > percent_adults, ]$reduction = NA}
  if (percent_twentyplus < 50) {df_Herzog[df_Herzog$strat == "twentyplus" & df_Herzog$vax_avail > percent_twentyplus, ]$reduction = NA}
  
  # 
  # p <- p + geom_line(data = df_Herzog[df_Herzog$variable == "var",],
  #                    aes(x = vax_avail, y = reduction, col = strat),
  #                    size = 1.1, alpha = 1)
  return(p)
} 

plot_serotesting_over_vax_avail = function(outcome, strat_to_plot, mycol){
  # plots comparison of testing v no testing
  if (strat_to_plot == "all"){plot_title <- "All Ages"
                              this_color = col_all
  } else if (strat_to_plot == "kids"){plot_title <- "Children (0-19)"
                                      this_color = col_kids
  } else if (strat_to_plot == "young adults"){plot_title <- "Young Adults (20-49)"
                                              this_color = col_youngadults
  } else if (strat_to_plot == "elderly"){plot_title <- "Elderly (60+)"
                                          this_color = col_elderly
  } else if (strat_to_plot == "adults"){plot_title <- "Adults (20+)"
                                        this_color = col_adults
  }
  
  total_cases <- rep(NA, 510)
  total_deaths <- rep(NA, 510)
  count <- 1
  for (i in list_all_sero_notest){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_all_sero_test){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_kids_sero_notest){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_kids_sero_test){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_adults_sero_notest){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_adults_sero_test){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_elderly_sero_notest){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_elderly_sero_test){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_twentyplus_sero_notest){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_twentyplus_sero_test){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  
  num_strategies <- 5
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies*2))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list*2), rep("kids", num_per_list*2), 
             rep("adults", num_per_list*2), rep("elderly", num_per_list*2), 
             rep("twentyplus", num_per_list*2))
  temp <-  c(rep("notest", num_per_list), rep("test", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  theme_set(theme_minimal(base_size = 26))
  
  percent_kids <- round(sum(age_demo[1:2])*100)
  percent_elderly <- round(sum(age_demo[7:9])*100)
  percent_adults <- round(sum(age_demo[3:5])*100)
  percent_twentyplus <- round(sum(age_demo[3:9])*100)
  
  if (outcome == "cases"){
    baseline_cases_notest <- compute_total_cases(list_all_sero_notest$`0`)
    baseline_cases_test <- compute_total_cases(list_all_sero_test$`0`)
    
    temp <- c(rep(baseline_cases_notest, num_per_list), rep(baseline_cases_test, num_per_list))
    baseline_cases <- c(rep(temp, num_strategies))
    
    reduction_in_cases <- (1-(total_cases/baseline_cases))*100
    
    df <- data.frame(vax_avail, strat, reduction_in_cases, variable)
    # for Belgium sero
    df[df$strat == "kids" & vax_avail > 23 & variable == "notest", ]$reduction_in_cases = NA
    df[df$strat == "kids" & vax_avail > 22 & variable == "test", ]$reduction_in_cases = NA
    df[df$strat == "elderly" & vax_avail > 26 & variable == "notest", ]$reduction_in_cases = NA
    df[df$strat == "elderly" & vax_avail > 25 & variable == "test", ]$reduction_in_cases = NA
    df[df$strat == "adults" & vax_avail > 38 & variable == "notest", ]$reduction_in_cases = NA
    df[df$strat == "adults" & vax_avail > 36 & variable == "test", ]$reduction_in_cases = NA
    
    # ## for NY sero
    # df[df$strat == "kids" & vax_avail > percent_kids & variable == "notest", ]$reduction_in_cases = NA
    # df[df$strat == "kids" & vax_avail > 18 & variable == "test", ]$reduction_in_cases = NA
    # df[df$strat == "elderly" & vax_avail > percent_elderly & variable == "notest", ]$reduction_in_cases = NA
    # df[df$strat == "elderly" & vax_avail > 18 & variable == "test", ]$reduction_in_cases = NA
    # df[df$strat == "adults" & vax_avail > percent_adults & variable == "notest", ]$reduction_in_cases = NA
    # df[df$strat == "adults" & vax_avail > 30 & variable == "test", ]$reduction_in_cases = NA
    # df[df$strat == "all" & vax_avail > 37 & variable == "test", ]$reduction_in_cases = NA
    
    #df <- df[df$strat == strat_to_plot, ]
    
    ggplot(df, aes(x = vax_avail, y = reduction_in_cases, col = strat)) +
      geom_line(aes(linetype = variable), size = 1.2, alpha = 0.9) +
      xlab("Total vaccine supply (% of pop)") +
      ylab("Reduction in infections (%)") +
      scale_linetype_discrete(name = "Scenario",
                              labels = c("w/o sero tests", "w/ sero tests"))+
      theme(
            #legend.text = element_text(size=22),
            #legend.title = element_text(size=22),
            legend.position = "none",
             axis.title.x = element_blank(),
             axis.text.y  = element_blank(),
             axis.text.x = element_text(size = 17),
             axis.title.y = element_blank(),
             plot.title = element_blank()) +
      scale_y_continuous(expand = c(0,0), limit = c(0, 100)) +
      scale_x_continuous(expand = c(0,0), limit = c(0, 50), breaks = c(0,10, 20, 30, 40, 50)) +
      scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                         labels =  c("Adults 20-49", "All Ages", "Adults 60+", 
                                     "Under 20", "Adults 20+", "Optimal")) +
      guides(colour = guide_legend(override.aes = list(size=6)))
  } 
  else if (outcome == "deaths"){
    baseline_deaths_notest <- compute_total_deaths(list_all_sero_notest$`0`)
    baseline_deaths_test <- compute_total_deaths(list_all_sero_test$`0`)
    
    temp <- c(rep(baseline_deaths_notest, num_per_list), rep(baseline_deaths_test, num_per_list))
    baseline_deaths <- c(rep(temp, num_strategies))
    
    reduction_in_deaths <- (1-(total_deaths/baseline_deaths))*100
    
    df <- data.frame(vax_avail, strat, reduction_in_deaths, variable)
    
    # for Belgium sero
    df[df$strat == "kids" & vax_avail > 23 & variable == "notest", ]$reduction_in_deaths = NA
    df[df$strat == "kids" & vax_avail > 22 & variable == "test", ]$reduction_in_deaths = NA
    df[df$strat == "elderly" & vax_avail > 26 & variable == "notest", ]$reduction_in_deaths = NA
    df[df$strat == "elderly" & vax_avail > 25 & variable == "test", ]$reduction_in_deaths = NA
    df[df$strat == "adults" & vax_avail > 38 & variable == "notest", ]$reduction_in_deaths = NA
    df[df$strat == "adults" & vax_avail > 36 & variable == "test", ]$reduction_in_deaths = NA
    
    # ## for NY sero
    # df[df$strat == "kids" & vax_avail > percent_kids & variable == "notest", ]$reduction_in_deaths = NA
    # df[df$strat == "kids" & vax_avail > 18 & variable == "test", ]$reduction_in_deaths = NA
    # df[df$strat == "elderly" & vax_avail > percent_elderly & variable == "notest", ]$reduction_in_deaths = NA
    # df[df$strat == "elderly" & vax_avail > 18 & variable == "test", ]$reduction_in_deaths = NA
    # df[df$strat == "adults" & vax_avail > percent_adults & variable == "notest", ]$reduction_in_deaths = NA
    # df[df$strat == "adults" & vax_avail > 30 & variable == "test", ]$reduction_in_deaths = NA
    # df[df$strat == "all" & vax_avail > 37 & variable == "test", ]$reduction_in_deaths = NA
    #df <- df[df$strat == strat_to_plot, ]
    
    ggplot(df, aes(x = vax_avail, y = reduction_in_deaths, col = strat)) +
      geom_line(aes(linetype = variable), size = 1.2, alpha = 0.9) +
      xlab("Total vaccine supply (% of pop)") +
      ylab("Reduction in deaths (%)") +
      scale_linetype_discrete(name = "Scenario",
                              labels = c("W/o sero tests", "W/ sero tests")) +
      scale_y_continuous(expand = c(0,0), limit = c(0, 100)) +
      scale_x_continuous(expand = c(0,0), limit = c(0, 50), breaks = c(0,10, 20, 30, 40, 50)) +
      theme(legend.position = "none", 
             axis.title.x = element_blank(),
             axis.text.x  = element_blank(),
             axis.text.y  = element_blank(),
             axis.title.y = element_blank(),
            plot.title = element_text(color = this_color, size = 17)) +
      scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                         labels =  c("Young Adults (20-49)", "All Ages", "Elderly (60+)", 
                                     "Children (0-19)", "Adults (20+)", "Optimal")) +
      guides(colour = guide_legend(override.aes = list(size=6))) 
      #ggtitle(plot_title)
  }
}

# Other plotting fn ----
get_legend = function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
