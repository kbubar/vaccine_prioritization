library(grid)
library(dplyr)
library(data.table)

move_vaccinated_event <- function(t, x, parms){
  v_e <- parms$v_e
  v_e_type <- parms$v_e_type
  num_perday <- parms$num_perday
  vax_proportion <- parms$vax_proportion
  groups <- parms$groups
  sp <- parms$sp
  se <- parms$se
  pop_total <- parms$pop_total
  
  # move those who are vaccinated in a given day
  num_compartment <- 13
  num_groups <- (length(x)-1)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  Sx   <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  E    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  Ev   <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Ex   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  I    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Iv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  Ix   <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  R    <- as.matrix(x[(9*num_groups+1):(10*num_groups)])
  Rv   <- as.matrix(x[(10*num_groups+1):(11*num_groups)])
  Rx   <- as.matrix(x[(11*num_groups+1):(12*num_groups)])
  D    <- as.matrix(x[(12*num_groups+1):(13*num_groups)])
  vax_supply <- x[13*num_groups+1]
  
  if (vax_supply >= num_perday*pop_total){
    nvax <- num_perday*pop_total
    vax_supply <- vax_supply - num_perday*pop_total
  } else {
    nvax <- vax_supply
    vax_supply <- 0
  }
  
  vax_distribution <- nvax*vax_proportion
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
    # don't go over the number eligible
    if (any(vax_distribution > vax_eligible)){
      vax_distribution[vax_distribution > vax_eligible] = vax_eligible
    }
  }
  
  alpha <- vax_distribution/vax_eligible
  alpha[alpha == Inf] <- 0 # no people in S,E,R
  alpha[is.nan(alpha)] <- 0 # no vax left and no people in S,E,R
  if(any(alpha > 1)){print("WARNING: alpha > 1 in move_vaccinated")
    alpha[alpha>1] <- 0} # more vaccines avaliable than eligible people
  
  if (v_e_type == "leaky") {
    dS  <- -(S*alpha*sp) - (S*alpha*(1-sp))
    dSv <- (S*alpha*sp)
    dSx <- (S*alpha*(1-sp))
    
    dE  <- -(E*alpha*sp) - (E*alpha*(1-sp))
    dEv <- (E*alpha*sp)
    dEx <- (E*alpha*(1-sp)) 
  } else {
    # all-or-nothing
    dS  <- -(S*alpha*sp) - (S*alpha*(1-sp))
    dSv <- (S*alpha*sp*v_e)
    dSx <- (S*alpha*sp*(1-v_e)) + (S*alpha*(1-sp)) 
    
    dE  <- - (E*alpha*sp) - (E*alpha*(1-sp))
    dEv <- (E*alpha*sp*v_e)
    dEx <- (E*alpha*sp*(1-v_e)) + (E*alpha*(1-sp))
  }
  
  dR <- - (R*alpha*(1-se)) - (R*alpha*se)
  dRv <- (R*alpha*(1-se))
  dRx <- (R*alpha*se)
  
  # update compartments
  S    <- S + dS
  Sv   <- Sv + dSv
  Sx   <- Sx + dSx
  E    <- E + dE
  Ev   <- Ev + dEv
  Ex   <- Ex + dEx
  R    <- R + dR
  Rv   <- Rv + dRv
  Rx   <- Rx + dRx
  
  # output updated compartments
  out <- c(S,Sv,Sx,E,Ev,Ex,I,Iv,Ix,R,Rv,Rx,D,vax_supply)
}

calculate_derivatives = function(t, x, parameters){
  # x is a vector of length (# model compartment types)*(# age groups)
  # S, E, I, R etc. are vectors of length num_groups
  num_compartment <- 13
  num_groups <- (length(x)-1)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  Sx   <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  E    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  Ev   <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Ex   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  I    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Iv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  Ix   <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  R    <- as.matrix(x[(9*num_groups+1):(10*num_groups)])
  Rv   <- as.matrix(x[(10*num_groups+1):(11*num_groups)])
  Rx   <- as.matrix(x[(11*num_groups+1):(12*num_groups)])
  D    <- as.matrix(x[(12*num_groups+1):(13*num_groups)])
  vax_supply <- x[13*num_groups+1]
  
  S[S < .Machine$double.eps] <- 0
  Sv[Sv < .Machine$double.eps] <- 0
  Sx[Sx < .Machine$double.eps] <- 0
  E[E < .Machine$double.eps] <- 0
  Ev[Ev < .Machine$double.eps] <- 0
  Ex[Ex < .Machine$double.eps] <- 0
  I[I < .Machine$double.eps] <- 0
  Iv[Iv < .Machine$double.eps] <- 0
  Ix[Ix < .Machine$double.eps] <- 0
  R[R < .Machine$double.eps] <- 0
  Rv[Rv < .Machine$double.eps] <- 0
  Rx[Rx < .Machine$double.eps] <- 0
  D[D < .Machine$double.eps] <- 0
  
  u <- parameters$u
  C <- parameters$C
  d_E <- parameters$d_E
  d_I <- parameters$d_I
  v_e <- parameters$v_e
  v_e_type <- parameters$v_e_type
  
  lambda <- C%*%((I+Iv+Ix)/(N_i-D))*u
  
  if (v_e_type == "leaky") {
    dSv <- -(Sv*(1-v_e)*lambda)
    dEv <- (Sv*(1-v_e)*lambda) - d_E*Ev 
  } else {
    # all-or-nothing
    dSv <- rep(0, num_groups)
    dEv <- -d_E*as.matrix(Ev)
  }
    
  dS  <- -(S*lambda)
  dSx <- -(Sx*lambda)

  dE  <- S*lambda - d_E*E
  dEx <- Sx*lambda - d_E*Ex

  dI  <- E*d_E - I*d_I
  dIv <- Ev*d_E - Iv*d_I
  dIx <- Ex*d_E - Ix*d_I

  dR  <- I*d_I*(1-IFR)
  dRv <- Iv*d_I*(1-IFR)
  dRx <- Ix*d_I*(1-IFR)

  dD  <- I*d_I*IFR + Iv*d_I*IFR + Ix*d_I*IFR
  
  dvax_supply <- 0
  
  out <- c(dS,dSv,dSx,dE,dEv,dEx,dI,dIv,dIx,dR,dRv,dRx,dD,dvax_supply)
  list(out)
}

calculate_derivatives_NTB = function(t, x, parameters){
  # x is a vector of length (# model compartment types)*(# age groups)
  # S, E, I, R etc. are vectors of length num_groups
  num_compartment <- 13
  num_groups <- (length(x)-1)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  Sx   <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  E    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  Ev   <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Ex   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  I    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Iv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  Ix   <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  R    <- as.matrix(x[(9*num_groups+1):(10*num_groups)])
  Rv   <- as.matrix(x[(10*num_groups+1):(11*num_groups)])
  Rx   <- as.matrix(x[(11*num_groups+1):(12*num_groups)])
  D    <- as.matrix(x[(12*num_groups+1):(13*num_groups)])
  vax_supply <- x[13*num_groups+1]
  
  S[S < .Machine$double.eps] <- 0
  Sv[Sv < .Machine$double.eps] <- 0
  Sx[Sx < .Machine$double.eps] <- 0
  E[E < .Machine$double.eps] <- 0
  Ev[Ev < .Machine$double.eps] <- 0
  Ex[Ex < .Machine$double.eps] <- 0
  I[I < .Machine$double.eps] <- 0
  Iv[Iv < .Machine$double.eps] <- 0
  Ix[Ix < .Machine$double.eps] <- 0
  R[R < .Machine$double.eps] <- 0
  Rv[Rv < .Machine$double.eps] <- 0
  Rx[Rx < .Machine$double.eps] <- 0
  D[D < .Machine$double.eps] <- 0
  
  u <- parameters$u
  C <- parameters$C
  ve_S <- parameters$ve_S
  ve_I <- parameters$ve_I
  ve_P <- parameters$ve_P
  
  lambda <- C %*% ((I+(1-ve_I)*Iv + Ix)/(N_i-D)) * u
  lambda_V <- (1-ve_S)*C %*% ((I+(1-ve_I)*Iv + Ix)/(N_i-D)) * u

  dSv <- -Sv*lambda_V
  dEv <- Sv*lambda_V - d_E*Ev
  
  dS  <- -S*lambda
  dE  <- S*lambda - d_E*E
  
  dSx <- -Sx*lambda
  dEx <- Sx*lambda - d_E*Ex
  
  dI  <- E*d_E - I*d_I
  dIv <- Ev*d_E - Iv*d_I
  dIx <- Ex*d_E - Ix*d_I
  
  dR  <- I*d_I*(1-IFR) 
  dRv <- Iv*d_I*(1 - (1-ve_P)*IFR)
  dRx <- Ix*d_I*(1-IFR)
  
  dD  <- I*d_I*IFR + Iv*d_I*(1-ve_P)*IFR + Ix*d_I*IFR
  
  dvax_supply <- 0
  
  out <- c(dS,dSv,dSx,dE,dEv,dEx,dI,dIv,dIx,dR,dRv,dRx,dD,dvax_supply)
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

get_best_strat_deaths_new = function(supply, all, kids, adults, 
                                     elderly, twentyplus){
  # supply is the vaccine supply of interest as an integer
  
  total_deaths <- c(compute_total_deaths_new(all[[supply + 1]]),
                    compute_total_deaths_new(kids[[supply + 1]]), 
                    compute_total_deaths_new(adults[[supply + 1]]),
                    compute_total_deaths_new(elderly[[supply + 1]]), 
                    compute_total_deaths_new(twentyplus[[supply + 1]]))
  
  baseline_deaths <- compute_total_deaths_new(all$`0`)
  reduction_in_deaths <- (1 - (total_deaths/baseline_deaths))*100
  
  strats <- c("all", "kids", "adults", "elderly", "twentyplus")
  strats[which.max(reduction_in_deaths)]
}

get_best_strat_cases_new = function(supply, all, kids, adults, 
                                    elderly, twentyplus){
  # supply is the vaccine supply of interest as an integer
  
  total_cases <- c(compute_total_cases_new(all[[supply + 1]]),
                   compute_total_cases_new(kids[[supply + 1]]), 
                   compute_total_cases_new(adults[[supply + 1]]),
                   compute_total_cases_new(elderly[[supply + 1]]), 
                   compute_total_cases_new(twentyplus[[supply + 1]]))
  
  baseline_cases <- compute_total_cases_new(all$`0`)
  reduction_in_cases <- (1 - (total_cases/baseline_cases))*100
  
  strats <- c("all", "kids", "adults", "elderly", "twentyplus")
  strats[which.max(reduction_in_cases)]
}

get_best_strat_yll_new = function(supply, all, kids, adults, 
                                  elderly, twentyplus){
  # supply is the vaccine supply of interest as an integer
  
  total_YLL <- c(compute_total_YLL_new(all[[supply + 1]]),
                 compute_total_YLL_new(kids[[supply + 1]]), 
                 compute_total_YLL_new(adults[[supply + 1]]),
                 compute_total_YLL_new(elderly[[supply + 1]]), 
                 compute_total_YLL_new(twentyplus[[supply + 1]]))
  
  baseline_YLL <- compute_total_YLL_new(all$`0`)
  reduction_in_YLL <- (1 - (total_YLL/baseline_YLL))*100
  
  strats <- c("all", "kids", "adults", "elderly", "twentyplus")
  strats[which.max(reduction_in_YLL)]
}

run_v_e_var = function(p, v_e_baseline, hinge_age, vax_amount, num_perday, v_e_type){
  this_v_e <- get_v_e(p, y0 = v_e_baseline, hinge_age = hinge_age) # hinge age is age group that v_e begins decreasing after
  
  all <- kids <- adults <- elderly <- twentyplus <- vector(mode = "list")
  
  if (vax_amount > 0){
    for (i in seq(0, vax_amount, by = vax_amount)){
      j <- vax_amount/100
      all[[paste0(i)]] <- run_sim_new(this_C, j, "all", num_perday, v_e_type, this_v_e)
      kids[[paste0(i)]] <- run_sim_new(this_C, j, "kids", num_perday, v_e_type, this_v_e)
      adults[[paste0(i)]] <- run_sim_new(this_C, j, "adults", num_perday, v_e_type, this_v_e)
      elderly[[paste0(i)]] <- run_sim_new(this_C, j, "elderly", num_perday, v_e_type, this_v_e)
      twentyplus[[paste0(i)]] <- run_sim_new(this_C, j, "twentyplus", num_perday, v_e_type, this_v_e)
    }
  } else{
      i <- 0
      all[[paste0(i)]] <- run_sim_new(this_C, vax_amount, "all", num_perday, v_e_type, this_v_e)
      kids[[paste0(i)]] <- run_sim_new(this_C, vax_amount, "kids", num_perday, v_e_type, this_v_e)
      adults[[paste0(i)]] <- run_sim_new(this_C, vax_amount, "adults", num_perday, v_e_type, this_v_e)
      elderly[[paste0(i)]] <- run_sim_new(this_C, vax_amount, "elderly", num_perday, v_e_type, this_v_e)
      twentyplus[[paste0(i)]] <- run_sim_new(this_C, vax_amount, "twentyplus", num_perday, v_e_type, this_v_e)
      i <- i + 1
      all[[paste0(i)]] <- all[[1]]
      kids[[paste0(i)]] <- kids[[1]]
      adults[[paste0(i)]] <- adults[[1]]
      elderly[[paste0(i)]] <- elderly[[1]]
      twentyplus[[paste0(i)]] <- twentyplus[[1]]
  }
  return(get_best_strat_for_tipping_point(all, kids, adults, elderly, twentyplus))
}

v_e_bisection = function(v_e_baseline, hinge_age, vax_amount, num_perday, v_e_type){
  p_high <- 1
  p_low <- 0
  
  temp <- run_v_e_var(p_high, v_e_baseline, hinge_age, vax_amount, num_perday, v_e_type)
  if (temp != "elderly"){
    print("Warning: elderly is not most strategic to start")
    return(list(100, temp))
  }
  
  temp <- run_v_e_var(p_low, v_e_baseline, hinge_age, vax_amount, num_perday, v_e_type)
  if (temp == "elderly"){
    print("Warning: elderly is always most strategic")
    return(list(0, temp))
  }
  
  running = TRUE
  while(running){
    p_next = (p_high + p_low)/2
    temp <- run_v_e_var(p_next, v_e_baseline, hinge_age, vax_amount, num_perday, v_e_type)
    if (temp == "elderly"){
      p_high <- p_next
    } else {
      p_low <- p_next
    }
    
    if ((p_high - p_low) < 0.01){
      strat <- run_v_e_var(p_low, v_e_baseline, hinge_age, vax_amount, num_perday, v_e_type)
      out <- list(p_next*100, strat)
      return(out)
    }
  }
}

get_best_strat_for_tipping_point = function(all, kids, adults, elderly, twentyplus){
  # supply is the vaccine supply of interest as an integer
  
  total_deaths <- c(compute_total_deaths_new(all[[2]]),
                    compute_total_deaths_new(kids[[2]]), 
                    compute_total_deaths_new(adults[[2]]),
                    compute_total_deaths_new(elderly[[2]]), 
                    compute_total_deaths_new(twentyplus[[2]]))
  
  baseline_deaths <- compute_total_deaths_new(all$`0`)
  reduction_in_deaths <- (1 - (total_deaths/baseline_deaths))*100
  
  strats <- c("all", "kids", "adults", "elderly", "twentyplus")
  strats[which.max(reduction_in_deaths)]
}

compute_R0 = function(u, C){
  gamma <- 1/5 # recovery period (I -> R), ref: Davies
  
  # Davies NGM
  Du <- diag(u, 9)
  Dy <- diag(1/gamma, 9)
  
  NGM <- Du %*% C %*% Dy
  R0  <- abs(eigen(NGM)$values[1])
}

compute_R_with_sero = function(u, C, sero){
  gamma <- 1/5 # recovery period (I -> R), ref: Davies
  
  # Davies NGM
  Du <- diag(u, 9)
  Dy <- diag(1/gamma, 9)
  D_sero <- diag(1 - sero)
  
  NGM <- D_sero %*% Du %*% C %*% Dy
  R  <- abs(eigen(NGM)$values[1])
}

scale_u_for_R0 = function(u, C, wanted_R0){
  scalehigh <- 10
  scalelow <- 100
  
  R0_high <- compute_R0(u/scalehigh, C)
  R0_low <- compute_R0(u/scalelow, C)
  
  if (R0_high < wanted_R0 || R0_low > wanted_R0){
    print("Error: wanted_R0 is not in the original range of R0 values")
    return(-1)
  }
  
  running = TRUE
  while(running){
    scale_next = (scalehigh + scalelow)/2
    temp <- compute_R0(u/scale_next, C)
    if (temp > wanted_R0){
      scalehigh <- scale_next
    } else {
      scalelow <- scale_next
    }
    
    if (abs(temp - wanted_R0) < 0.0001){
      return(scale_next)
    }
  }
}

compute_total_cases_new = function(df){
  infections <- rep(0,num_groups) # number of age groups
  
  R_index <- 83 # col number for R
  Rv_index <- 92 # col number for Rv
  Rx_index <- 101 # col number for Rx
  D_index <- 110
  
  for (i in 1:num_groups) {
    # infections = total # recovered - initial recovered (seropositive)
    infections[i] <- df[dim(df)[1], R_index] - df[1, R_index] +
      df[dim(df)[1], Rv_index] - df[1, Rv_index] + 
      df[dim(df)[1], Rx_index] - df[1, Rx_index] + 
      df[dim(df)[1], D_index]
    R_index  <- R_index + 1
    Rv_index <- Rv_index + 1
    Rx_index <- Rx_index + 1
    D_index  <- D_index + 1
  }
  tot_infections <- sum(infections)/pop_total * 100
}

compute_total_deaths_new = function(df){
  deaths <- rep(0,num_groups)
  D_index <- 110
  
  for (i in 1:num_groups) {
    deaths[i] <- df[dim(df)[1], D_index]
    D_index <- D_index + 1
  }
  
  tot_deaths <- sum(deaths)/pop_total * 100
}

compute_total_YLL_new = function(df){
  YLL <- rep(0,num_groups)
  
  D_index <- 110
  for (i in 1:num_groups) {
    deaths_temp <- df[dim(df)[1], D_index]
    YLL[i] <- YLL_vec[i]*deaths_temp
    D_index <- D_index + 1
  }
  
  tot_YLL <- sum(YLL) # absolute number of YLL
}

compile_total_cases = function(total){
  count <- 1
  for (i in list_all){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  
  total
}

compile_total_cases_var = function(total){
  count <- 1
  for (i in list_all){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  
  total
}

compile_total_deaths = function(total){
  count <- 1
  for (i in list_all){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  
  total
}

compile_total_deaths_var = function(total){
  count <- 1
  for (i in list_all){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  total
}

compile_total_YLL = function(total){
  count <- 1
  for (i in list_all){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  
  total
}

compile_total_YLL_var = function(total){
  count <- 1
  for (i in list_all){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  total
}

get_reduction_df = function(outcome){
  total <- rep(NA, 255)
  
  if (outcome == "cases"){
    total <- compile_total_cases(total)
    baseline <- compute_total_cases_new(list_all$`0`)
  } else if (outcome == "deaths"){
    total <- compile_total_deaths(total)
    baseline <- compute_total_deaths_new(list_all$`0`)
  } else if (outcome == "YLL"){
    total <- compile_total_YLL(total)
    baseline <- compute_total_YLL_new(list_all$`0`)
  }

  num_strategies <- 5
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list), rep("kids", num_per_list), 
             rep("adults", num_per_list), rep("elderly", num_per_list), 
             rep("twentyplus", num_per_list))
  variable <-  c(rep("constant", num_per_list))
  
  
  baseline <- c(rep(baseline, num_per_list))
  
  reduction <- (1-(total/baseline))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

get_reduction_df_var = function(outcome){
  total <- rep(NA, 510)
  
  if (outcome == "cases"){
    total <- compile_total_cases_var(total)
    baseline <- compute_total_cases_new(list_all$`0`)
    baseline_var <- compute_total_cases_new(list_all_var$`0`)
  } else if (outcome == "deaths"){
    total <- compile_total_deaths_var(total)
    baseline <- compute_total_deaths_new(list_all$`0`)
    baseline_var <- compute_total_deaths_new(list_all_var$`0`)
  } else if (outcome == "YLL"){
    total <- compile_total_YLL_var(total)
    baseline <- compute_total_YLL_new(list_all$`0`)
    baseline_var <- compute_total_YLL_new(list_all_var$`0`)
  }
  
  num_strategies <- 5
  vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies*2))
  num_per_list <- 51
  strat <- c(rep("all", num_per_list*2), 
             rep("kids", num_per_list*2), 
             rep("adults", num_per_list*2),
             rep("elderly", num_per_list*2), 
             rep("twentyplus", num_per_list*2))
  temp <-  c(rep("constant", num_per_list), rep("var", num_per_list))
  variable <- c(rep(temp, num_strategies))
  
  temp <- c(rep(baseline, num_per_list), rep(baseline_var, num_per_list))
  baseline <- c(rep(temp, num_strategies))
  
  reduction <- (1-(total/baseline))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

# Plot trajectories ----
get_rowindex <- function(df, rName) {
  which(match(colnames(df), rName) == 1)
  }

gather_compartments_overtime <- function(df, compartment, strat){
  final_time <- as.numeric(dim(df)[1])-1
  
  if (compartment == "I"){
    start <- "I1"
    end <- "Ix9"
  } else if (compartment == "R"){
    start <- "R1"
    end <- "D9"
  } else if (compartment == "D"){
    start <- "D1"
    end <- "D9"
  }
  
  total_df <- data.frame(time = 0:final_time,
                         percent = ((apply(df[, get_rowindex(df, start):get_rowindex(df, end)], 
                                        1, sum))/pop_total)*100,
                         strat = rep(strat, final_time + 1))
}

plot_strat_overtime = function(compartment, df_baseline, df_all, df_adults, df_kids, df_twentyplus, 
                               df_elderly, this_time) {
  df_baseline <- gather_compartments_overtime(df_baseline, compartment, "no vax")
  df_all <- gather_compartments_overtime(df_all, compartment, "all")
  df_adults <- gather_compartments_overtime(df_adults, compartment, "adults")
  df_kids <- gather_compartments_overtime(df_kids, compartment, "kids")
  df_twentyplus <- gather_compartments_overtime(df_twentyplus, compartment, "twentyplus")
  df_elderly <- gather_compartments_overtime(df_elderly, compartment, "elderly")
  
  df <- rbind(df_adults, df_all, df_elderly, df_kids, df_twentyplus)
  
  p <- ggplot(df, aes(x = time, y = percent)) +
    theme_minimal(base_size = 12) +
    annotate("rect", xmin=0, xmax=this_time, ymin=0, ymax=Inf, alpha=0.5, fill= "#dddddd") +
    geom_line(data = df_baseline, aes(x = time, y = percent), col = "#4F4F4F", size = 0.5, 
              linetype = "dashed") +
    geom_line(aes(color = strat), size = 0.5) +
    xlab("Time (days)") +
    scale_x_continuous(expand = c(0,0), limit = c(0, 250)) +#, breaks = c(0, 100, 200, 300)) +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Adults 20-49", "All Ages", "Adults 60+", 
                                   "Under 20", "Adults 20+")) +
    theme(legend.position = "none") 
  
  if (compartment == "I") {
    ymax <- 4
    p <- p + ylab("\nInfected (%)")
  } else if (compartment == "R") {
    ymax <- 60
    p <- p + ylab("Cumulative\nincidence (%)")
  } else if (compartment == "D") {
    ymax <- 0.6
    p <- p + ylab("Cumulative\nmortality (%)") 
  }
  
  p <- p + scale_y_continuous(expand = c(0,0), limit = c(0, ymax)) + 
    coord_fixed(250*4/(5*ymax))
}

gather_compartments_byage <- function(df, compartment, age_group){
  final_time <- as.numeric(dim(df)[1])-1
  
  if (compartment == "I"){
    this_comp   <- paste0("I", age_group)
    this_comp_v <- paste0("Iv", age_group)
    this_comp_x <- paste0("Ix", age_group)
    
    total_df <- data.frame(time = 0:final_time,
                           percent = ((df[, get_rowindex(df, this_comp)] + 
                                               df[, get_rowindex(df, this_comp_v)] +
                                               df[, get_rowindex(df, this_comp_x)])/pop_total)*100,
                           age_group = rep(age_group, final_time + 1))
  } else if (compartment == "R"){
    this_comp   <- paste0("R", age_group)
    this_comp_v <- paste0("Rv", age_group)
    this_comp_x <- paste0("Rx", age_group)
    this_comp_d <- paste0("D", age_group)
    
    total_df <- data.frame(time = 0:final_time,
                           percent = ((df[, get_rowindex(df, this_comp)] + 
                                               df[, get_rowindex(df, this_comp_v)] +
                                               df[, get_rowindex(df, this_comp_x)] + 
                                               df[, get_rowindex(df, this_comp_d)])/pop_total)*100,
                           age_group = rep(age_group, final_time + 1))
  } else if (compartment == "D"){
    this_comp <- paste0("D", age_group)
    
    total_df <- data.frame(time = 0:final_time,
                           percent = ((df[, get_rowindex(df, this_comp)])/pop_total)*100,
                           age_group = rep(age_group, final_time + 1))
  }
}

plot_allages_onestrategy = function(this_df, compartment, this_time, groups, this_col){
  # INPUTS:
  # df: simulation to plot
  # compartment: character of compartment to plot i.e. "I", "R"
  #
  # OUTPUT:
  # p: ggplot object of plot
  
  df_toplot <- data.frame()
  for (i in 1:9) {
    temp <- gather_compartments_byage(this_df, compartment, i)
    df_toplot <- rbind(df_toplot, temp)
  }
  
  # Plot
  p <- ggplot(df_toplot[!(df_toplot$age_group %in% groups),], aes(x = time, y = percent)) +
    annotate("rect", xmin=0, xmax=this_time, ymin=0, ymax=Inf, alpha=0.5, fill="#dddddd") +
    geom_line(aes(group = factor(age_group)), size = 0.4, col = "#cccccc") +
    xlab("Time (Days)") +
    scale_x_continuous(expand = c(0,0),limit = c(0, 250), breaks = c(0,100,200)) +
    # col_grey + 
    theme(legend.position = "none") 
  
  p <- p + geom_line(data = df_toplot[df_toplot$age_group %in% groups,],
                     aes(x = time, y = percent, group = factor(age_group)),
                     color = this_col, size = 0.4)
  if (compartment == "I") {
    ymax <- 0.7
    p <- p + ylab("\nInfected (%)") + 
      scale_y_continuous(expand = c(0,0), limit = c(0, ymax))#E, breaks = c(0, 1, 2))
  } else if (compartment == "R") {
    ymax <- 10
    p <- p + ylab("Cumulative\nincidence (%)") + 
      scale_y_continuous(expand = c(0,0), limit = c(0, ymax), breaks = c(0, 5, 10, 15))
  } else if (compartment == "D") {
    ymax <- 0.3
    p <- p + ylab("Cumulative\nmortality (%)") + 
      scale_y_continuous(expand = c(0,0), limit = c(0, ymax))#, breaks = c(0, 0.25, 0.5))
  }
  p <- p + coord_fixed(250*2/(3.5*ymax))
}

plot_supp_dynamics_panel = function(compartment){
  t_one <- 11
  infect_10 <- plot_strat_overtime(compartment, list_all[[1]], list_all[[t_one]], list_adults[[t_one]], 
                                   list_kids[[t_one]], list_twentyplus[[t_one]], list_elderly[[t_one]], 50) +
    alllabels_theme +
    theme(axis.title.x = element_blank()) + 
    theme(axis.title.y = element_blank()) + 
    ggtitle("10% vaccine supply") + 
    theme(plot.title = element_text(color = "black"))
  
  t_two <- 31
  infect_30 <- plot_strat_overtime(compartment, list_all[[1]], list_all[[t_two]], list_adults[[t_two]], 
                                   list_kids[[t_two]], list_twentyplus[[t_two]], list_elderly[[t_two]], 150) + 
    onlyx_theme + 
    theme(axis.title.x = element_blank()) + 
    ggtitle("30% vaccine supply") + 
    theme(plot.title = element_text(color = "black"))
  
  ###
  shading1 <- 50
  shading2 <- 150
  infect_none1 <- plot_allages_onestrategy(list_elderly[[1]], compartment, 0, NA, "#cccccc") + 
    onlyy_theme +
    theme(axis.title.y = element_blank())
  infect_adults1 <- plot_allages_onestrategy(list_adults[[t_one]], compartment, shading1, 3:5, col_youngadults) + 
    onlyy_theme + 
    theme(axis.title.y = element_blank())
  infect_elderly1 <- plot_allages_onestrategy(list_elderly[[t_one]], compartment, shading1, 7:9, col_elderly) + 
    alllabels_theme +
    theme(axis.title.x = element_blank()) + 
    theme(axis.title.y = element_blank())
  infect_all1 <- plot_allages_onestrategy(list_all[[t_one]], compartment, shading1, 1:9, col_all) + 
    onlyx_theme +
    theme(axis.title.x = element_blank()) 
  infect_twentyplus1 <- plot_allages_onestrategy(list_twentyplus[[t_one]], compartment, shading1, 3:9, col_adults) + 
    nolabels_theme
  infect_kids1 <- plot_allages_onestrategy(list_kids[[t_one]], compartment, shading1, 1:2, col_kids) + 
    nolabels_theme
  
  infect_none2 <- plot_allages_onestrategy(list_elderly[[1]], compartment, 0, NA, "#cccccc") + 
    nolabels_theme
  infect_adults2 <- plot_allages_onestrategy(list_adults[[t_two]], compartment, shading2, 3:5, col_youngadults) + 
    nolabels_theme
  infect_elderly2 <- plot_allages_onestrategy(list_elderly[[t_two]], compartment, shading2, 7:9, col_elderly) + 
    alllabels_theme +
    theme(axis.title.x = element_blank()) + 
    theme(axis.title.y = element_blank()) + 
    theme(axis.text.y = element_blank())
  infect_all2 <- plot_allages_onestrategy(list_all[[t_two]], compartment, shading2, 1:9, col_all) + 
    onlyx_theme +
    theme(axis.title.x = element_blank()) 
  infect_twentyplus2 <- plot_allages_onestrategy(list_twentyplus[[t_two]], compartment, shading2, 3:9, col_adults) + 
    nolabels_theme
  infect_kids2 <- plot_allages_onestrategy(list_kids[[t_two]], compartment, shading2, 1:2, col_kids) + 
    nolabels_theme
  
  dynamics <- ggarrange(infect_none1, infect_kids1, infect_none2, infect_kids2,
                             infect_adults1, infect_twentyplus1, infect_adults2, infect_twentyplus2,
                             infect_elderly1, infect_all1, infect_elderly2, infect_all2,
                             nrow = 3)

  panel_infect_strat <- ggarrange(infect_10, infect_30,
                                  nrow = 1)
  # export as 8x4.25
  if (compartment == "I") {
    grid.arrange(panel_infect_strat, dynamics,
                 heights = c(1, 1),
                 bottom = "Time (days)",
                 left = "Infected (% of total population)",
                 right = "")
  } else if (compartment == "R"){
    grid.arrange(panel_infect_strat, dynamics,
                 heights = c(1, 1),
                 bottom = "Time (days)",
                 left = "Cumulative incidence (% of total population)",
                 right = "")
  } else {
    grid.arrange(panel_infect_strat, dynamics,
                 heights = c(1, 1),
                 bottom = "Time (days)",
                 left = "Cumulative mortality since initial vaccine rollout\n(% of total population)",
                 right = "")
  }
}

# Bar plots ----
plot_vaxdist_hist = function(){
  p1 <- barplot_vax_strat("kids") + 
    theme(axis.title.y = element_blank())
  p2 <- barplot_vax_strat("adults") + 
    theme(axis.title.y = element_blank())
  p3 <- barplot_vax_strat("20+") + 
    theme(axis.title.y = element_blank())
  p4 <- barplot_vax_strat("elderly") + 
    theme(axis.title.y = element_blank())
  p5 <- barplot_vax_strat("all") + 
    theme(axis.title.y = element_blank())
  
  panel <- ggarrange(p1, p2, p3, p4, p5,
                     nrow = 5,
                     left = textGrob("Distribution of vaccines (%)", rot = 90, hjust = 0.5),
                     bottom = textGrob("Age (years)", vjust = 0))
}

barplot_vax_strat = function(strategy){
  if (strategy == "no vax"){
    vaccinated <- rep(0,num_groups) * 100
    this_color <- "#808080"
    plot_title <- "No Vaccines"
  } else {
    if (strategy == "all"){
      groups <- 1:9
      this_color <- col_all
      plot_title <- "All ages"
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
  #age_groups <- c("0","10", "20", "30", "40", "50", "60", "70", "80")
  age_groups <- seq(0, 80, by = 10)
  vax_proportion <- vax_proportion*100
  
  df <- data.frame(age_groups, vax_proportion, age_demo)
  
  # plot
  p <- ggplot(df, aes(x=age_groups)) + 
    #geom_bar(aes(y=age_demo), fill="white",color="black", position="stack", stat="identity") + 
    geom_bar(aes(y=vax_proportion), position="stack", stat="identity", fill = this_color) + 
    xlab("Age group") +
    ggtitle(plot_title) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 52), breaks = c(0, 25, 50), labels = c(0, "", 50))
  
  if (strategy == "all"){
    p <- p + scale_x_continuous(limits = c(-5, NA), expand = c(0,0),
                                breaks=seq(-5, 75, by = 10),
                                labels=c("0", "10", "20", "30", "40", "50", "60", "70", "80")) +
      theme(plot.title = element_text(color = this_color, size = 10),#, face = "bold"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size = 9))
  } else {
    p <- p + theme(plot.title = element_text(color = this_color,  size = 10),
                   axis.title.x =element_blank(),
                   axis.text.x = element_blank(),
                   axis.title.y = element_blank()) +
      scale_x_continuous(limits = c(-5, NA), expand = c(0,0),
                         breaks=seq(-5, 75, by = 10),
                         labels=c("", "", "", "", "", "", "", "", ""))
  }
  p <- p + theme(panel.grid.minor = element_blank(),
                 panel.grid.major = element_blank(),
                 axis.ticks = element_line(size = 0.35),
                 axis.line = element_line(size = 0.35),
                 plot.title = element_text(vjust = -1),
                 plot.margin = unit(c(-0.06,0,0,0), "cm"))
  return(p)
}

when_strat_switch = function(list_df, groups){
  # find the vaccine supply where you switch from strategy to vaccinating everyone else
  # input: list_df is a list with df for a strategy for vax supply 0-50
  #        groups: vector index of the strategy
  # output: vax supply where switch occurs
  
  other_groups <- 1:9
  other_groups <- other_groups[!other_groups %in% groups]
  x_switch <- -1
  
  for (i in 1:length(list_df)){
    temp <- list_df[[i]]
    if (any(temp[dim(temp)[1], other_groups + 10] > 0)){
      x_switch <- i-1 
      break
    }
  }
  return(x_switch)
}

plot_over_vax_avail = function(outcome, var = FALSE){
  library(ggplot2)
  theme_set(theme_minimal(base_size = 12))
  
  x_adults_switch     <- when_strat_switch(list_adults, 3:5)
  x_kids_switch       <- when_strat_switch(list_kids, 1:2)
  x_elderly_switch    <- when_strat_switch(list_elderly, 7:9)
  x_twentyplus_switch <- when_strat_switch(list_twentyplus, 3:9)
  
  # get dataframe for specific outcome
  if (var){
    df <- get_reduction_df_var(outcome)
    
    x_adults_switch_var     <- when_strat_switch(list_adults_var, 3:5)
    x_kids_switch_var       <- when_strat_switch(list_kids_var, 1:2)
    x_elderly_switch_var    <- when_strat_switch(list_elderly_var, 7:9)
    x_twentyplus_switch_var <- when_strat_switch(list_twentyplus_var, 3:9)
  } else {df <- get_reduction_df(outcome)}
  
  points_x <- c(x_kids_switch, x_elderly_switch)
  points_y <- c(df[df$strat == "kids" & df$vax_avail == x_kids_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "elderly" & df$vax_avail == x_elderly_switch & df$variable == "constant", ]$reduction)
  
  if (x_adults_switch > 0){
    points_x <- c(points_x, x_adults_switch)
    points_y <- c(points_y, df[df$strat == "adults" & df$vax_avail == x_adults_switch & df$variable == "constant", ]$reduction)
  }
  if (x_twentyplus_switch > 0){
    points_x <- c(points_x, x_twentyplus_switch)
    points_y <- c(points_y, df[df$strat == "twentyplus" & df$vax_avail == x_twentyplus_switch & df$variable == "constant", ]$reduction)
  }
  
  if (var){
    points_x <- c(points_x, x_kids_switch_var, x_elderly_switch_var)
    points_y <- c(points_y, df[df$strat == "kids" & df$vax_avail == x_kids_switch_var & df$variable == "var", ]$reduction,
                  df[df$strat == "elderly" & df$vax_avail == x_elderly_switch_var & df$variable == "var", ]$reduction)
    
    if (x_adults_switch_var > 0){
      points_x <- c(points_x, x_adults_switch_var)
      points_y <- c(points_y, df[df$strat == "adults" & df$vax_avail == x_adults_switch_var & df$variable == "var", ]$reduction)
    }
    if (x_twentyplus_switch_var > 0){
      points_x <- c(points_x, x_twentyplus_switch_var)
      points_y <- c(points_y, df[df$strat == "twentyplus" & df$vax_avail == x_twentyplus_switch_var & df$variable == "var", ]$reduction)
    }
  }
  
  df_var <- df[df$variable == "var", ]
  df_const <- df[df$variable == "constant", ]

  # plot
  p <- ggplot() +
    geom_line(aes(x = df_const$vax_avail, y = df_const$reduction, col = df_const$strat),
              size = 0.5, alpha = 0.9) +
    geom_line(aes(x = df_var$vax_avail, y = df_var$reduction, col = df_var$strat), 
              linetype = "dashed", size = 1, alpha = 0.9) +
    xlab("Total vaccine supply (% of population)") +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Adults 20-49", "All Ages", "Adults 60+", 
                                   "Under 20", "Adults 20+")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 75), breaks = c(0, 25, 50, 75)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    coord_fixed(50*4/(5*75)) +
    theme(legend.position = "none") +
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  p <- p + geom_point(aes(x = points_x, y = points_y), size = 1) 
  
  if (outcome == "cases"){ p <- p + ylab("Reduction in\ninfections (%)")}
  else if (outcome == "deaths"){ p <- p + ylab("Reduction in\ndeaths (%)")}
  else if (outcome == "YLL"){ p <- p + ylab("Reduction in\nYLL (%)")}
  return(p)
}

# Other plotting fn ----
plot_age_dep_ve = function() {
# age-dependent ve plot
  groups <- c(0,10,20,30,40,50,60,70,80,90)
  v_e_plot <- v_e_var
  v_e_plot[10] <- v_e_var[9]
  
  df <- data.frame(groups, v_e_plot)
  
  age_dep_ve <- ggplot(df, aes(x = groups, y = v_e_plot*100)) +
    geom_step(size = 0.8, linetype = "dashed") +
    geom_hline(yintercept = 90, size = 0.5) +
    geom_point(aes(x = 60, y = 90), size = 2) +
    geom_point(aes(x = 89.5, y = 50), size = 2) +
    ylab("Efficacy (%)") +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 90.2), breaks = c(0,20,40,60,80),
                       labels = c(0,20,40,60,80)) +
    #theme(panel.grid.minor = element_blank()) +
    xlab("Age (years)") +
    theme(#axis.text = element_text(size = 9),
          axis.title.y = element_text(vjust = 0.2)) +
    coord_fixed(90*4/(5*100)) 
}

plot_heatmap = function(data, num, param){
  
  p <- ggplot(data, aes(x = factor(Var2), y = factor(Var1), fill = factor(value)))+
    geom_tile() + 
    scale_fill_manual(breaks = seq(0, 5, by = 1),  
                      values = c("#E6E6E6", col_youngadults, col_adults, col_elderly, col_kids, col_all)) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text = element_text(size = 10)) +
    scale_x_discrete(expand = c(0,0), breaks = seq(0, 50, by = 10)) + 
    coord_fixed(50*4/(6*num))
  
  if (param == "R0"){
    p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 2),
                                labels = seq(1.1, 2.0, by = 0.1))
  } else if (param == "ve"){
    p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 2),
                                labels = c(" 30", " 50", " 70", " 90"))
  } else if (param == "country") {
    p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 1),
                              labels = c("BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"))
  } else if (param == "rollout") {
    p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 1),
                              labels = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1))
  } else if (param == "ve_I"){
    p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 2),
                              labels = c(" 0"," 20", " 40", " 60", " 80", "100"))
  }
  
  
  for (i in seq(1.5, num + 0.5, by = 1)){
    p <- p + geom_hline(yintercept = i, col = "white")
  }

  for (i in seq(10.5, 40.5, by = 10)){
    p <- p + geom_vline(xintercept = i, col = "white")
  }
  
  p
}

plot_tipping_point_heatmap = function(list_tp){
  baselines <- c(rep(baseline_vec[1], 50),
                 rep(baseline_vec[2], 50),
                 rep(baseline_vec[3], 50),
                 rep(baseline_vec[4], 50),
                 rep(baseline_vec[5], 50),
                 rep(baseline_vec[6], 50))
  vaxsupply <- c(rep(1:50, 6))
  
  strats <- c(list_tp[[1]][[2]],
              list_tp[[2]][[2]],
              list_tp[[3]][[2]],
              list_tp[[4]][[2]],
              list_tp[[5]][[2]],
              list_tp[[6]][[2]])
  
  df <- data.frame(vaxsupply, baselines, strats)
  
  p <- ggplot(df, aes(x = vaxsupply, y = baselines*100, fill = factor(strats)))+
    geom_tile() + 
    # scale_fill_manual(breaks = seq(1, 5, by = 1),  
    #                   values = c(col_youngadults, col_adults, col_elderly, col_kids, col_all)) +
    colFill + 
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text = element_text(size = 10)) +
    #axis.title.y = element_blank()) +
    scale_x_continuous(expand = c(0,0), breaks = seq(0, 50, by = 10)) +
    scale_y_continuous(expand = c(0,0), breaks = seq(50, 100, by = 10)) 

  p
}

plot_tipping_point_heatmap_2 = function(list_tp){
  baselines <- c(rep(baseline_vec[1], 50),
                 rep(baseline_vec[2], 50),
                 rep(baseline_vec[3], 50),
                 rep(baseline_vec[4], 50),
                 rep(baseline_vec[5], 50),
                 rep(baseline_vec[6], 50))
  vaxsupply <- c(rep(1:50, 6))
  
  strats <- c(list_tp[[1]][[1]],
              list_tp[[2]][[1]],
              list_tp[[3]][[1]],
              list_tp[[4]][[1]],
              list_tp[[5]][[1]],
              list_tp[[6]][[1]])
  
  breakpoints <- c(-1, 0.000001, 25, 50, 99.999999, 101)
  names <- c("None", "0-25%", "25-50%", "50-100%", "NA")
  
  df <- data.frame(vaxsupply, baselines, strats)
  df$cat <- cut(df$strats, breaks = breakpoints, labels = names)
  
  p <- ggplot(df, aes(x = factor(vaxsupply), y = factor(baselines*100), fill = cat))+
    geom_tile() + 
    colFill2 +
    theme(axis.title.x = element_blank()) +
          #axis.text = element_text(size = 10)) +
    scale_x_discrete(expand = c(0,0), breaks = seq(0, 50, by = 10)) +
    scale_y_discrete(expand = c(0,0), breaks = seq(50, 100, by = 10))+
    coord_fixed(50*4/(6*6))
  
  for (i in seq(1.5, 5 + 0.5, by = 1)){
    p <- p + geom_hline(yintercept = i, col = "white")
  }
  
  for (i in seq(10.5, 40.5, by = 10)){
    p <- p + geom_vline(xintercept = i, col = "white")
  }
  
  p
}

plot_opt_histogram = function(this_df, this_percent_vax){
  
  dist_of_vax <- matrix(nrow = 21, ncol = 9)
  
  for (i in 1:21){
    tot <- sum(N_i*this_df[i, 2:10])
    dist_of_vax[i,] <- as.numeric(((N_i*this_df[i, 2:10])/tot)*100)
  }
  
  dist_of_vax <- as.data.frame(dist_of_vax)
  colnames(dist_of_vax) <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  
  dist_of_vax$Percent_vax <- 0:20
  
  new_df <- dist_of_vax %>%
    select(Percent_vax, "0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+") %>%
    gather(key = "age_group", value = "num", -Percent_vax)
  
  p <- ggplot(new_df[new_df$Percent_vax == this_percent_vax,], aes(y=num, x=age_group, 
                                                                   fill = as.factor(Percent_vax), 
                                                                   group = as.factor(Percent_vax))) + 
    geom_bar(position = "stack", stat = "identity", fill = "#E6AB02") +
    #xlab("Age group") + 
    #ylab("Age distribution of vaccines (%)")+
    theme_classic(base_size = 12) +
    scale_x_discrete(breaks=c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"), 
                     labels=c("", "", "", "", "", "", "", "", "")) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 75), breaks = c(0,25,50,75)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x =element_blank(),
          axis.title.y =element_blank())
  
  return(p)
}

get_tp_values = function(list_tp){
  vals <- matrix(NA, nrow = 6, ncol = 50)
  vals[1,] <- list_tp[[1]][[1]]
  vals[2,] <- list_tp[[2]][[1]]
  vals[3,] <- list_tp[[3]][[1]]
  vals[4,] <- list_tp[[4]][[1]]
  vals[5,] <- list_tp[[5]][[1]]
  vals[6,] <- list_tp[[6]][[1]]
  vals
}

get_legend = function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
