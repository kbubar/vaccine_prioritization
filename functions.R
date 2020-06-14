library(grid)

calculate_derivatives=function(t, x, parameters){
  # the parameters in the parameters list are:
  #    the probability of transmission on contact, beta
  #    the incubation period, nu
  #    the recovery period, gamma
  #    the contact matrix, C, that is the # contacts per day among age groups
  #
  # Note that x is a vector of length (#model compartment types)*(#age classes)
  # Thus, S, E, I and R are vectors, all of length nage
  ncompartment <- 5
  nage <- length(x)/ncompartment
  S    <- as.matrix(x[1:nage])
  E    <- as.matrix(x[(nage+1):(2*nage)])
  I    <- as.matrix(x[(2*nage+1):(3*nage)])
  R    <- as.matrix(x[(3*nage+1):(4*nage)])
  V    <- as.matrix(x[(4*nage+1):(5*nage)])
  
  I[I<0] = 0
  with(as.list(parameters),{
    # note that because S, I and R are all vectors of length nage, so will N,
    # and dS, dI, and dR
    beta = c(beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9)
    
    N = S+E+I+R+V
    dS = -as.matrix(S*beta)*(as.matrix(C)%*%as.matrix(I/N))
    dE = -dS - nu*as.matrix(E)
    dI = nu*as.matrix(E) - gamma*as.matrix(I)
    dR = +gamma*as.matrix(I)
    dV = 0*as.matrix(I)
    
    out=c(dS,dE,dI,dR,dV)
    list(out)
  })
}

compute_total_cases = function(df){
  infections <- rep(0,nage)
  
  R_index <- 29 # col number for R1
  
  for (i in 1:nage) {
    infections[i] <- max(df[R_index])
    R_index <- R_index + 1
  }
  
  tot_infections <- sum(infections)/npop * 100
  tot_infections <- round(tot_infections, 1)
}

compute_total_deaths = function(df){
  deaths <- rep(0,nage)
  
  R_index <- 29 # col number for R1
  
  for (i in 1:nage) {
    deaths[i] <- max(df[R_index])*IFR[i]
    R_index <- R_index + 1
  }
  
  tot_deaths <- sum(deaths)/npop * 100
  tot_deaths <- round(tot_deaths, 3)
}

# Plot trajectories ----
plot_allages_onestrategy =function(df, compartment){
  # INPUTS:
  # df: simulation to plot
  # compartment: character of compartment to plot i.e. "I", "R"
  #
  # OUTPUT:
  # p: ggplot object of plot
  
  # make new df for plotting
  col_to_gather <- vector(mode="character", length=nage)
  for (i in 1:9) {col_to_gather[i] <- paste0(compartment,i)}
  
  new_df <- df %>%
    select(time, col_to_gather) %>%
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
  
  # Plot
  theme_set(theme_minimal(base_size = 15))
  
  string_to_add <- paste0("Total infected: ", compute_total_cases(df), "%")
  
  # Create a text
  grob <- grobTree(textGrob(string_to_add, x=0.7,  y=0.95, hjust=0,
                            gp=gpar(col="red", fontsize=13, fontface="italic")))
  # Plot
  
  
  p <- ggplot(new_df, aes(x = time, y = percent)) + 
    geom_line(aes(color = age_group), size = 2) +
    xlab("Time") +
    annotation_custom(grob) +
    scale_color_brewer(palette = "Spectral", name = "Age Group", 
                       labels =  c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"))
  
  if (compartment == "I") {
    p <- p + ylab("Percent Infected") + ylim(0,0.4)
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
plot_finalT_onestrategy = function(solution_df, IFR, vaccinated){
  infections <- rep(0,nage)
  deaths <- rep(0, nage)
  final_S <- rep(0, nage)
  
  R_index <- 29 # col number for R1

  for (i in 1:nage) {
    deaths[i] <- max(solution_df[R_index])*IFR[i]/N[i]
    infections[i] <- max(solution_df[R_index])*(1-IFR[i])/N[i] # recovered that didn't die
    final_S[i] <- solution_df[(max(t) + 1), (i + 1)]/N[i]
    R_index <- R_index + 1
  }
  
  age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  
  df <- data.frame(age_groups, final_S, deaths, infections, vaccinated)
  
  new_df <- df %>%
    # select(age_groups, infections, deaths, vaccinated, final_S) %>%
    select(age_groups, vaccinated) %>%
    gather(key = "type", value = "num", -age_groups)
  
  ggplot(new_df, aes(fill=type, y=num, x=age_groups)) + 
    geom_bar(position="stack", stat="identity", fill = col3) + 
    ylab("% Vaccinated") + 
    xlab("Age group") + 
    ylim(0,1)
  
  # ggplot(df, aes(x=age_groups, y=deaths)) + 
  #   geom_bar(stat = "identity") + 
  #   ylim(0, 0.1)
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

# Other plotting fn ----
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
