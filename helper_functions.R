library(grid)

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
 
  u <- parameters$beta
  C <- parameters$C
  nu <- parameters$nu
  gamma <- parameters$gamma
  v_e <- parameters$v_e
  
  N = S+E+I+R+V
  dS = -as.matrix(S*u)*(as.matrix(C)%*%as.matrix(I/N))
  dV = (1-v_e)*(-as.matrix(V*u))*(as.matrix(C)%*%as.matrix(I/N))
  dE = -dS - dV - nu*as.matrix(E)
  dI = nu*as.matrix(E) - gamma*as.matrix(I)
  dR = +gamma*as.matrix(I)
  

  
  out=c(dS,dE,dI,dR,dV)
  list(out)
}


get_v_e = function(p){
  # INPUT: p is the final v_e % (as a decimal) 
  #       i.e. p = 1 -> perfect vaccine at all ages
  # ASSUMPTIONS: perfect v_e up to age 50 then linear decline to p
  # OUTPUT: vector v_e 
  
  x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
  y <- ((p-1)/30)*(x - 50) + 1
  y[1:6] <- 1
  
  y
}

compute_R0 = function(u){
  # Davies NGM
  Du <- diag(u, 9)
  Dy <- diag(1/gamma, 9)
  
  NGM <- Du %*% C %*% Dy
  R0  <- abs(eigen(NGM)$values[1])
}

compute_total_cases = function(df){
  infections <- rep(0,num_groups)
  
  R_index <- 29 # col number for R1
  
  for (i in 1:num_groups) {
    infections[i] <- max(df[R_index])
    R_index <- R_index + 1
  }
  
  tot_infections <- sum(infections)/pop_total * 100
  tot_infections <- round(tot_infections, 1)
}

compute_total_deaths = function(df){
  deaths <- rep(0,num_groups)
  
  R_index <- 29 # col number for R1
  
  for (i in 1:num_groups) {
    deaths[i] <- max(df[R_index])*IFR[i]
    R_index <- R_index + 1
  }
  
  tot_deaths <- sum(deaths)/pop_total * 100
  #tot_deaths <- round(tot_deaths, 3)
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
  
  # Plot
  theme_set(theme_minimal(base_size = 15))
  
  string_to_add <- paste0("Total infected: ", compute_total_cases(df), "%")
  
  # Create a text
  grob <- grobTree(textGrob(string_to_add, x=0.5,  y=0.95, hjust=0,
                            gp=gpar(col="white", fontsize=13, fontface="bold")))
  # grob <- grobTree(textGrob(string_to_add, just= c("right", "top"),
  #                           gp=gpar(col="white", fontsize=13, fontface="bold")))
  
  # Plot
  
  
  p <- ggplot(new_df, aes(x = time, y = percent)) + 
    theme_dark() + 
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
barplot_at_finalT = function(solution_df, IFR, strategy, y){
  if (strategy == "no vax"){
    vaccinated <- rep(0,num_groups) * 100
    this_color <- "#808080"
    plot_title <- "No Vaccines"
  } 
  else if (strategy == "all"){
    vaccinated <- vax_propall/N_i * 100
    this_color <- col_all
    plot_title <- "Vaccinate All"
  } 
  else if (strategy == "kids"){
    vaccinated <- vax_propkids/N_i * 100
    this_color <- col_kids
    plot_title <- "Vaccinate Kids"
  } 
  else if (strategy == "adults"){
    vaccinated <- vax_propadults/N_i * 100
    this_color <- col_adults
    plot_title <- "Vaccinate Adults"
  } 
  else if (strategy == "elderly"){
    vaccinated <- vax_propelderly/N_i * 100
    vaccinated[vaccinated>100] = 100
    this_color <- col_elderly
    plot_title <- "Vaccinate Elderly"
  }
  
  infections <- rep(0,num_groups)
  deaths <- rep(0, num_groups)
  final_S <- rep(0, num_groups)
  
  R_index <- 29 # col number for R1

  for (i in 1:num_groups) {
    deaths[i] <- max(solution_df[R_index])*IFR[i]/N_i[i]*100
    infections[i] <- max(solution_df[R_index])/N_i[i]*100 # total infected (including dead)
    # infections[i] <- max(solution_df[R_index])*(1-IFR[i])/N_i[i] # total infected (excluding dead)
    final_S[i] <- solution_df[(max(t) + 1), (i + 1)]/N_i[i]
    R_index <- R_index + 1
  }
  
  age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  
  df <- data.frame(age_groups, final_S, deaths, infections, vaccinated)
  
  new_df <- df %>%
    # select(age_groups, infections, deaths, vaccinated, final_S) %>%
    select(age_groups, all_of(y)) %>%
    gather(key = "type", value = "num", -age_groups)
  
  theme_set(theme_minimal(base_size = 12))
  
  # plot
  p <- ggplot(new_df, aes(fill=type, y=num, x=age_groups)) + 
    geom_bar(position="stack", stat="identity", fill = this_color) + 
    xlab("Age group") + 
    scale_x_discrete(breaks=c("0-9","20-29", "40-49", "60-69", "80+")) + 
    theme(plot.subtitle = element_text(size = 9, face = "italic", hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  if (y == "infections") {
      string_to_add <- paste0("Total infected: ", compute_total_cases(solution_df), "%")

      p <- p + ylim(0, 100) + 
        ylab("Infected (%)") + 
        # labs(title = plot_title, subtitle = string_to_add)
        labs(subtitle = string_to_add)
    } 
  else if (y == "deaths") {
      string_to_add <- paste0("Total deaths: ", round(compute_total_deaths(solution_df),2), "%")
      
      p <- p + ylim(0, 8) + 
        ylab("Deaths (%)") + 
        # labs(title = plot_title, subtitle = string_to_add)
        labs(subtitle = string_to_add)
    } 
  else if (y == "vaccinated"){
      p <- p + ylim(0, 100) + 
        ylab("Vaccinated (%)") + 
        labs(title = plot_title)
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
plot_over_vax_avail = function(outcome, var_name, df_novax_var, list_all_var, list_kids_var, list_adults_var, list_elderly_var){
  total_cases <- rep(NA, 408)
  total_deaths <- rep(NA, 408)
  
  count <- 1
  for (i in list_all){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total_cases[count] <- compute_total_cases(i)
    total_deaths[count] <- compute_total_deaths(i)
    count <- count + 1
  }
  
  # # To 50 by 5s
  #vax_avail <- c(rep(seq(0, 50, by = 5), 8))
  # strat <- c(rep("all", 22), rep("kids", 22), rep("adults", 22), rep("elderly", 22))
  # num_per_list <- 11
  # susceptibility <- c(rep("constant", num_per_list), rep("var", num_per_list), rep("constant", num_per_list), rep("var", num_per_list), rep("constant", num_per_list), rep("var", num_per_list), rep("constant", num_per_list), rep("var", num_per_list))
  
  # # to 50 by 1s
  vax_avail <- c(rep(seq(0, 50, by = 1), 8))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list*2), rep("kids", num_per_list*2), rep("adults", num_per_list*2), rep("elderly", num_per_list*2))
  variable <- c(rep("constant", num_per_list), rep("var", num_per_list), rep("constant", num_per_list), rep("var", num_per_list), rep("constant", num_per_list), rep("var", num_per_list), rep("constant", num_per_list), rep("var", num_per_list))
  
  baseline_cases <- compute_total_cases(df_novax)
  baseline_deaths <- compute_total_deaths(df_novax)
  baseline_cases_var <- compute_total_cases(df_novax_var)
  baseline_deaths_var <- compute_total_deaths(df_novax_var)
  
  baseline_cases <- c(rep(baseline_cases, num_per_list), rep(baseline_cases_var, num_per_list), rep(baseline_cases, num_per_list), rep(baseline_cases_var, num_per_list), rep(baseline_cases, num_per_list), rep(baseline_cases_var, num_per_list), rep(baseline_cases, num_per_list), rep(baseline_cases_var, num_per_list))
  baseline_deaths <- c(rep(baseline_deaths, num_per_list), rep(baseline_deaths_var, num_per_list), rep(baseline_deaths, num_per_list), rep(baseline_deaths_var, num_per_list), rep(baseline_deaths, num_per_list), rep(baseline_deaths_var, num_per_list), rep(baseline_deaths, num_per_list), rep(baseline_deaths_var, num_per_list))
  
  reduction_in_cases <- (1-(total_cases/baseline_cases))*100
  reduction_in_deaths <- (1 - (total_deaths/baseline_deaths))*100
  
  df <- data.frame(vax_avail, strat, reduction_in_cases, total_deaths, reduction_in_deaths, variable)
  
  if (outcome == "cases"){
    ggplot(df, aes(x = vax_avail, y = reduction_in_cases, col = strat, fill = strat)) +
      geom_abline(slope = 1, intercept = 0, size = 3, alpha = 0.2) +
      geom_line(aes(linetype = variable), size = 2, alpha = 0.9) +
      xlab("Vaccine available (% of total pop)") +
      ylab("% Reduction in cases") +
      #ggtitle("Optimizing for least cases") +
      scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                         labels =  c("Adults", "All", "Elderly", "Kids")) +
      scale_fill_brewer(palette = "Dark2", name = "Allocation Strategy",
                         labels =  c("Adults", "All", "Elderly", "Kids")) +
      scale_linetype_discrete(name = var_name,
                              labels = c("Constant", "Age-dependent"))+
      ylim(0, 100) +
      theme(legend.position = "none")
  } 
  else if (outcome == "deaths"){
    ggplot(df, aes(x = vax_avail, y = reduction_in_deaths, col = strat, fill = strat)) +
        geom_abline(slope = 1, intercept = 0, size = 3, alpha = 0.2) +
        geom_line(aes(linetype = variable), size = 2, alpha = 0.9) +
        xlab("Vaccine available (% of total pop)") +
        ylab("% Reduction in deaths") +
        #ggtitle("Optimizing for least deaths") +
        scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                           labels =  c("Adults", "All", "Elderly", "Kids")) +
        scale_fill_brewer(palette = "Dark2", name = "Allocation Strategy",
                          labels =  c("Adults", "All", "Elderly", "Kids")) +
        scale_linetype_discrete(name = var_name,
                                labels = c("Constant", "Age-dependent")) +
        ylim(0, 100) 
        # theme(legend.position = "none")
  }
  reduction_in_cases
}

# Other plotting fn ----
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
