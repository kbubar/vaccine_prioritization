library(grid)
library(dplyr)
library(data.table)
# calculate_derivatives=function(t, x, parameters){
#   # the parameters in the parameters list are:
#   #    the probability of transmission on contact, beta
#   #    the incubation period, nu
#   #    the recovery period, gamma
#   #    the contact matrix, C, that is the # contacts per day among age groups
#   #
#   # Note that x is a vector of length (#model compartment types)*(#age classes)
#   # Thus, S, E, I and R are vectors, all of length num_groups
#   ncompartment <- 5
#   num_groups <- length(x)/ncompartment
#   S    <- as.matrix(x[1:num_groups])
#   E    <- as.matrix(x[(num_groups+1):(2*num_groups)])
#   I    <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
#   R    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
#   V    <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
#   
#   I[I<0] = 0
#   
#   u <- parameters$u
#   C <- parameters$C
#   d_E <- parameters$d_E
#   d_I <- parameters$d_I
#   v_e <- parameters$v_e
#   
#   N = S+E+I+R+V
#   dS = -as.matrix(S*u)*(as.matrix(C)%*%as.matrix(I/N))
#   dV = 0*(1-v_e)*(-as.matrix(V*u))*(as.matrix(C)%*%as.matrix(I/N)) # multiply by 0 if all-or-nothing
#   dE = -dS - dV - d_E*as.matrix(E)
#   dI = d_E*as.matrix(E) - d_I*as.matrix(I)
#   dR = +d_I*as.matrix(I)
#   
#   out=c(dS,dE,dI,dR,dV)
#   list(out)
# }

### new framework ###
move_vaccinated = function(x, num_perday, vax_supply, vax_proportion, groups, v_e, 
                           v_e_type, sp, se) {
  # move those who are vaccinated in a given day
  num_compartment <- 13
  num_groups <- length(x)/num_compartment
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
  
  if (vax_supply >= num_perday*pop_total){
    nvax <- num_perday*pop_total
  } else {
    nvax <- vax_supply
  }
  
  vax_distribution <- nvax*vax_proportion
  vax_eligible <- (S+E)*sp + R*(1-se)
  #vax_distribution[vax_distribution > (S+E+R)] <- S[vax_distribution > (S+E+R)] + E[vax_distribution > (S+E+R)] + R[vax_distribution > (S+E+R)]
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
    # don't go over the number eligible
    if (any(vax_distribution > vax_eligible)){
      vax_distribution[vax_distribution > vax_eligible] = vax_eligible
      }
  }
  
  alpha <- vax_distribution/vax_eligible
  alpha[alpha == Inf] <- 0 # no people in S,E,R
  alpha[is.nan(alpha)] <- 0 # no vax left and no people in S,E,R
  #print(paste0("t: ", t, " max alpha:", max(alpha))) ### FIXME
  if(any(alpha > 1)){print("ERROR: alpha > 1 in move_vaccinated")}
  
  if (v_e_type == "leaky") {
    dS  <- -as.matrix(S*alpha*sp) - as.matrix(S*alpha*(1-sp))
    dSv <- as.matrix(S*alpha*sp)
    dSx <- as.matrix(S*alpha*(1-sp))
    
    dE  <- -as.matrix(E*alpha*sp) - as.matrix(E*alpha*(1-sp))
    dEv <- as.matrix(E*alpha*sp)
    dEx <- as.matrix(E*alpha*(1-sp)) 
  } else {
    # all-or-nothing
    dS  <- -as.matrix(S*alpha*sp) - as.matrix(S*alpha*(1-sp))
    dSv <- as.matrix(S*alpha*sp*v_e)
    dSx <- as.matrix(S*alpha*sp*(1-v_e)) + as.matrix(S*alpha*(1-sp)) 
    
    dE  <- - as.matrix(E*alpha*sp) - as.matrix(E*alpha*(1-sp))
    dEv <- as.matrix(E*alpha*sp*v_e)
    dEx <- as.matrix(E*alpha*sp*(1-v_e)) + as.matrix(E*alpha*(1-sp))
  }
  
  dR <- - as.matrix(R*alpha*(1-se)) - as.matrix(R*alpha*se)
  dRv <- as.matrix(R*alpha*(1-se))
  dRx <- as.matrix(R*alpha*se)
  
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
  out <- c(S,Sv,Sx,E,Ev,Ex,I,Iv,Ix,R,Rv,Rx,D)
}

move_vaccinated_NTB = function(x, num_perday, vax_supply, vax_proportion, groups) {
  # move those who are vaccinated in a given day
  num_compartment <- 9
  num_groups <- length(x)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  E    <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  Ev   <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  I    <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Iv   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  R    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Rv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  D    <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  
  if (vax_supply >= num_perday*pop_total){
    nvax <- num_perday*pop_total
  } else {
    nvax <- vax_supply
  }
  
  vax_distribution <- nvax*vax_proportion
  vax_eligible <- S+E+R
  #vax_distribution[vax_distribution > (S+E+R)] <- S[vax_distribution > (S+E+R)] + E[vax_distribution > (S+E+R)] + R[vax_distribution > (S+E+R)]
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
    # don't go over the number eligible
    if (any(vax_distribution > vax_eligible)){
      vax_distribution[vax_distribution > vax_eligible] = vax_eligible
    }
  }
  
  alpha <- vax_distribution/vax_eligible
  alpha[alpha == Inf] <- 0 # no people in S,E,R
  alpha[is.nan(alpha)] <- 0 # no vax left and no people in S,E,R
  #print(paste0("t: ", t, " max alpha:", max(alpha))) ### FIXME
  if(any(alpha > 1)){print("ERROR: alpha > 1 in move_vaccinated")}

  dS  <- -as.matrix(S*alpha)
  dSv <- as.matrix(S*alpha)
  
  dE  <- -as.matrix(E*alpha)
  dEv <- as.matrix(E*alpha)
  
  dR <- -as.matrix(R*alpha)
  dRv <- as.matrix(R*alpha)
  
  # update compartments
  S    <- S + dS
  Sv   <- Sv + dSv
  E    <- E + dE
  Ev   <- Ev + dEv
  R    <- R + dR
  Rv   <- Rv + dRv
  
  # output updated compartments
  out <- c(S,Sv,E,Ev,I,Iv,R,Rv,D)
}

calculate_derivatives_new=function(t, x, parameters){
  # x is a vector of length (# model compartment types)*(# age groups)
  # S, E, I, R etc. are vectors of length num_groups
  num_compartment <- 13
  num_groups <- length(x)/num_compartment
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
  
  # I[I<0] = 0
  # Iv[Iv<0] = 0
  # Ix[Ix<0] = 0
  
  u <- parameters$u
  C <- parameters$C
  d_E <- parameters$d_E
  d_I <- parameters$d_I
  v_e <- parameters$v_e
  v_e_type <- parameters$v_e_type
  
  lambda <- as.matrix((C)%*%as.matrix((I+Iv+Ix)/N_i) * as.matrix(u) )
  
  if (v_e_type == "leaky") {
    dSv <- -as.matrix(Sv*(1-v_e)*lambda)
    dEv <- as.matrix(Sv*(1-v_e)*lambda) - d_E*as.matrix(Ev) 
  } else {
    # all-or-nothing
    dSv <- rep(0, num_groups)
    dEv <- -d_E*as.matrix(Ev)
  }
  
  dS  <- -as.matrix(S*lambda)
  dSx <- -as.matrix(Sx*lambda)
  
  dE  <- as.matrix(S*lambda) - d_E*as.matrix(E)
  dEx <- as.matrix(Sx*lambda) - d_E*as.matrix(Ex)
  
  dI  <- as.matrix(E*d_E) - as.matrix(I*d_I)
  dIv <- as.matrix(Ev*d_E) - as.matrix(Iv*d_I)
  dIx <- as.matrix(Ex*d_E) - as.matrix(Ix*d_I)
  
  dR  <- as.matrix(I*d_I*(1-IFR)) 
  dRv <- as.matrix(Iv*d_I*(1-IFR)) 
  dRx <- as.matrix(Ix*d_I*(1-IFR)) 
  
  dD  <- as.matrix(I*d_I*IFR + Iv*d_I*IFR + Ix*d_I*IFR)
    
  out <- c(dS,dSv,dSx,dE,dEv,dEx,dI,dIv,dIx,dR,dRv,dRx,dD)
  list(out)
}

calculate_derivatives_NTB = function(t, x, parameters){
  # x is a vector of length (# model compartment types)*(# age groups)
  # S, E, I, R etc. are vectors of length num_groups
  num_compartment <- 9
  num_groups <- length(x)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  E    <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  Ev   <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  I    <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  Iv   <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  R    <- as.matrix(x[(6*num_groups+1):(7*num_groups)])
  Rv   <- as.matrix(x[(7*num_groups+1):(8*num_groups)])
  D    <- as.matrix(x[(8*num_groups+1):(9*num_groups)])
  
  S[S < .Machine$double.eps] <- 0
  Sv[Sv < .Machine$double.eps] <- 0
  E[E < .Machine$double.eps] <- 0
  Ev[Ev < .Machine$double.eps] <- 0
  I[I < .Machine$double.eps] <- 0
  Iv[Iv < .Machine$double.eps] <- 0
  R[R < .Machine$double.eps] <- 0
  Rv[Rv < .Machine$double.eps] <- 0
  D[D < .Machine$double.eps] <- 0
  
  u <- parameters$u
  C <- parameters$C
  ve_S <- parameters$ve_S
  ve_I <- parameters$ve_I
  ve_P <- parameters$ve_P
  
  lambda <- as.matrix(C %*% as.matrix(((I+(1-ve_I)*Iv)/N_i)) * as.matrix(u))
  lambda_V <- (1-ve_S)*as.matrix(C %*% as.matrix((I+(1-ve_I)*Iv)/N_i) * as.matrix(u))

  dSv <- -as.matrix(Sv*lambda_V)
  dEv <- as.matrix(Sv*lambda_V) - d_E*as.matrix(Ev) 
  
  dS  <- -as.matrix(S*lambda)
  dE  <- as.matrix(S*lambda) - d_E*as.matrix(E)
  
  dI  <- as.matrix(E*d_E) - as.matrix(I*d_I)
  dIv <- as.matrix(Ev*d_E) - as.matrix(Iv*d_I)
  
  dR  <- as.matrix(I*d_I*(1-IFR)) 
  dRv <- as.matrix(Iv*d_I*(1 - (1-ve_P)*IFR))
  
  dD  <- as.matrix(I*d_I*IFR + Iv*d_I*(1-ve_P)*IFR)
  
  out <- c(dS,dSv,dE,dEv,dI,dIv,dR,dRv,dD)
  list(out)
}

# calculate_derivatives_vax_overtime=function(t, x, parameters){
#   
#   ncompartment <- 5
#   num_groups <- length(x)/ncompartment
#   S    <- as.matrix(x[1:num_groups])
#   E    <- as.matrix(x[(num_groups+1):(2*num_groups)])
#   I    <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
#   R    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
#   V    <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
#   
#   I[I<0] = 0
#   
#   u <- parameters$u
#   C <- parameters$C
#   d_E <- parameters$d_E
#   d_I <- parameters$d_I
#   v_e <- parameters$v_e
#   vax_proportion <- parameters$vax_proportion
#   vax_supply <- parameters$vax_supply
#   npop <- parameters$npop
#   N <- parameters$N
#   num_perday <- parameters$num_perday
#   
#   if (vax_supply >= num_perday*npop){
#     nvax <- num_perday*npop
#   }
#   else {
#     nvax <- vax_supply
#   }
#   
#   vax_distribution <- nvax*vax_proportion
#   vax_distribution[vax_distribution > N] <- N[vax_distribution > N]
#   
#   prob_vaccinated <- vax_distribution/(S+E+R)
#   
#   # # leaky
#   # V = V + as.matrix(vax_distribution) # update V so that those who are vaccinating have ve% effectiveness on the day they're vaccinated
#   # dS = -as.matrix(S*u)*(as.matrix(C)%*%as.matrix(I/N)) - as.matrix(S*prob_vaccinated)
#   # dV = as.matrix(vax_distribution) - as.matrix(R*prob_vaccinated) - (as.matrix((1-v_e)*V*u))*(as.matrix(C)%*%as.matrix(I/N))
#   # dE = as.matrix(S*u)*(as.matrix(C)%*%as.matrix(I/N)) - d_E*as.matrix(E) - as.matrix(E*prob_vaccinated) + (1-v_e)*(as.matrix(V*u))*(as.matrix(C)%*%as.matrix(I/N)) 
#   # dI = d_E*as.matrix(E) - d_I*as.matrix(I)
#   # dR = d_I*as.matrix(I)
#   
#   # all-or-nothing
#   dS = -as.matrix(S*u)*(as.matrix(C)%*%as.matrix(I/N)) - as.matrix(S*prob_vaccinated*v_e)
#   dV = as.matrix(vax_distribution*v_e) - as.matrix(R*prob_vaccinated*v_e)
#   dE = as.matrix(S*u)*(as.matrix(C)%*%as.matrix(I/N)) - d_E*as.matrix(E) - as.matrix(E*prob_vaccinated*v_e)
#   dI = d_E*as.matrix(E) - d_I*as.matrix(I)
#   dR = d_I*as.matrix(I)
#   
#   out=c(dS,dE,dI,dR,dV)
#   list(out)
# }

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

# get_best_strat_deaths = function(supply, all, kids, adults, 
#                                  elderly, twentyplus){
#   # supply is the vaccine supply of interest as an integer
#   
#   total_deaths <- c(compute_total_deaths(all[[supply + 1]]),
#                     compute_total_deaths(kids[[supply + 1]]), 
#                     compute_total_deaths(adults[[supply + 1]]),
#                     compute_total_deaths(elderly[[supply + 1]]), 
#                     compute_total_deaths(twentyplus[[supply + 1]]))
#   
#   baseline_deaths <- compute_total_deaths(all$`0`)
#   reduction_in_deaths <- (1 - (total_deaths/baseline_deaths))*100
#   
#   strats <- c("all", "kids", "adults", "elderly", "twentyplus")
#   strats[which.max(reduction_in_deaths)]
# }
# 
# get_best_strat_cases = function(supply, all, kids, adults, 
#                                 elderly, twentyplus){
#   # supply is the vaccine supply of interest as an integer
#   
#   total_cases <- c(compute_total_cases(all[[supply + 1]]),
#                    compute_total_cases(kids[[supply + 1]]), 
#                    compute_total_cases(adults[[supply + 1]]),
#                    compute_total_cases(elderly[[supply + 1]]), 
#                    compute_total_cases(twentyplus[[supply + 1]]))
#   
#   baseline_cases <- compute_total_cases(all$`0`)
#   reduction_in_cases <- (1 - (total_cases/baseline_cases))*100
#   
#   strats <- c("all", "kids", "adults", "elderly", "twentyplus")
#   strats[which.max(reduction_in_cases)]
# }
# 
# get_best_strat_yll = function(supply, all, kids, adults, 
#                               elderly, twentyplus){
#   # supply is the vaccine supply of interest as an integer
#   
#   total_YLL <- c(compute_total_YLL(all[[supply + 1]]),
#                  compute_total_YLL(kids[[supply + 1]]), 
#                  compute_total_YLL(adults[[supply + 1]]),
#                  compute_total_YLL(elderly[[supply + 1]]), 
#                  compute_total_YLL(twentyplus[[supply + 1]]))
#   
#   baseline_YLL <- compute_total_YLL(all$`0`)
#   reduction_in_YLL <- (1 - (total_YLL/baseline_YLL))*100
#   
#   strats <- c("all", "kids", "adults", "elderly", "twentyplus")
#   strats[which.max(reduction_in_YLL)]
# }

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


get_best_strat_deaths_NTB = function(supply, all, kids, adults, 
                                     elderly, twentyplus){
  # supply is the vaccine supply of interest as an integer
  
  total_deaths <- c(compute_total_deaths_NTB(all[[supply + 1]]),
                    compute_total_deaths_NTB(kids[[supply + 1]]), 
                    compute_total_deaths_NTB(adults[[supply + 1]]),
                    compute_total_deaths_NTB(elderly[[supply + 1]]), 
                    compute_total_deaths_NTB(twentyplus[[supply + 1]]))
  
  baseline_deaths <- compute_total_deaths_NTB(all$`0`)
  reduction_in_deaths <- (1 - (total_deaths/baseline_deaths))*100
  
  strats <- c("all", "kids", "adults", "elderly", "twentyplus")
  strats[which.max(reduction_in_deaths)]
}

get_best_strat_cases_NTB = function(supply, all, kids, adults, 
                                    elderly, twentyplus){
  # supply is the vaccine supply of interest as an integer
  
  total_cases <- c(compute_total_cases_NTB(all[[supply + 1]]),
                   compute_total_cases_NTB(kids[[supply + 1]]), 
                   compute_total_cases_NTB(adults[[supply + 1]]),
                   compute_total_cases_NTB(elderly[[supply + 1]]), 
                   compute_total_cases_NTB(twentyplus[[supply + 1]]))
  
  baseline_cases <- compute_total_cases_NTB(all$`0`)
  reduction_in_cases <- (1 - (total_cases/baseline_cases))*100
  
  strats <- c("all", "kids", "adults", "elderly", "twentyplus")
  strats[which.max(reduction_in_cases)]
}

get_best_strat_yll_NTB = function(supply, all, kids, adults, 
                                  elderly, twentyplus){
  # supply is the vaccine supply of interest as an integer
  
  total_YLL <- c(compute_total_YLL_NTB(all[[supply + 1]]),
                 compute_total_YLL_NTB(kids[[supply + 1]]), 
                 compute_total_YLL_NTB(adults[[supply + 1]]),
                 compute_total_YLL_NTB(elderly[[supply + 1]]), 
                 compute_total_YLL_NTB(twentyplus[[supply + 1]]))
  
  baseline_YLL <- compute_total_YLL_NTB(all$`0`)
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
      all[[paste0(i)]] <- run_sim_new(C, j, "all", num_perday, v_e_type, this_v_e)
      kids[[paste0(i)]] <- run_sim_new(C, j, "kids", num_perday, v_e_type, this_v_e)
      adults[[paste0(i)]] <- run_sim_new(C, j, "adults", num_perday, v_e_type, this_v_e)
      elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", num_perday, v_e_type, this_v_e)
      twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", num_perday, v_e_type, this_v_e)
    }
  } else{
      i <- 0
      all[[paste0(i)]] <- run_sim_new(C, vax_amount, "all", num_perday, v_e_type, this_v_e)
      kids[[paste0(i)]] <- run_sim_new(C, vax_amount, "kids", num_perday, v_e_type, this_v_e)
      adults[[paste0(i)]] <- run_sim_new(C, vax_amount, "adults", num_perday, v_e_type, this_v_e)
      elderly[[paste0(i)]] <- run_sim_new(C, vax_amount, "elderly", num_perday, v_e_type, this_v_e)
      twentyplus[[paste0(i)]] <- run_sim_new(C, vax_amount, "twentyplus", num_perday, v_e_type, this_v_e)
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
    
    if ((p_high - p_low) < 0.001){
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
  scalehigh <- 30
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

# compute_total_cases = function(df){
#   infections <- rep(0,num_groups) # number of age groups
#   
#   R_index <- 29 # col number for R1
#   
#   for (i in 1:num_groups) {
#     # infections = total # recovered - initial recovered (seropositive)
#     infections[i] <- max(df[R_index]) - min(df[R_index])
#     R_index <- R_index + 1
#   }
#   
#   tot_infections <- sum(infections)/pop_total * 100
# }

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

compute_total_cases_NTB = function(df){
  infections <- rep(0,num_groups) # number of age groups
  
  R_index <- 56 # col number for R
  Rv_index <- 65 # col number for Rv
  D_index <- 74
  
  for (i in 1:num_groups) {
    # infections = total # recovered - initial recovered (seropositive)
    infections[i] <- df[dim(df)[1], R_index] - df[1, R_index] +
      df[dim(df)[1], Rv_index] - df[1, Rv_index] + 
      df[dim(df)[1], D_index]
    R_index  <- R_index + 1
    Rv_index <- Rv_index + 1
    D_index  <- D_index + 1
  }
  tot_infections <- sum(infections)/pop_total * 100
}
# compute_total_deaths = function(df){
#   deaths <- rep(0,num_groups)
#   
#   R_index <- 29 # col number for R1
#   
#   for (i in 1:num_groups) {
#     # infections = total # recovered - initial recovered (seropositive)
#     deaths[i] <- (max(df[R_index]) - min(df[R_index]))*IFR[i]
#     R_index <- R_index + 1
#   }
#   
#   tot_deaths <- sum(deaths)/pop_total * 100
# }

compute_total_deaths_new = function(df){
  deaths <- rep(0,num_groups)
  D_index <- 110
  
  for (i in 1:num_groups) {
    deaths[i] <- df[dim(df)[1], D_index]
    D_index <- D_index + 1
  }
  
  tot_deaths <- sum(deaths)/pop_total * 100
}

compute_total_deaths_NTB = function(df){
  deaths <- rep(0,num_groups)
  D_index <- 74
  
  for (i in 1:num_groups) {
    deaths[i] <- df[dim(df)[1], D_index]
    D_index <- D_index + 1
  }
  
  tot_deaths <- sum(deaths)/pop_total * 100
}

# compute_total_YLL = function(df){
#   YLL <- rep(0,num_groups)
#   
#   R_index <- 29 # col number for R1
#   
#   for (i in 1:num_groups) {
#     # infections = total # recovered - initial recovered (seropositive)
#     deaths_temp <- (max(df[R_index]) - min(df[R_index]))*IFR[i]
#     YLL[i] <- YLL_vec[i]*deaths_temp
#     R_index <- R_index + 1
#   }
#   
#   tot_YLL <- sum(YLL) # absolute number of YLL
# }

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

compute_total_YLL_NTB = function(df){
  YLL <- rep(0,num_groups)
  
  D_index <- 74
  for (i in 1:num_groups) {
    deaths_temp <- df[dim(df)[1], D_index]
    YLL[i] <- YLL_vec[i]*deaths_temp
    D_index <- D_index + 1
  }
  
  tot_YLL <- sum(YLL) # absolute number of YLL
}

# get_reduction_in_cases_df = function(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
#   total_cases <- rep(NA, 510)
#   
#   count <- 1
#   for (i in list_all){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_all_var){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_kids){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_kids_var){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_adults){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_adults_var){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_elderly){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_elderly_var){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_twentyplus){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   for (i in list_twentyplus_var){
#     total_cases[count] <- compute_total_cases(i)
#     count <- count + 1
#   }
#   
#   num_strategies <- 5
#   vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies*2))
#   num_per_list <- 51
#   strat <- c(rep("all", num_per_list*2), rep("kids", num_per_list*2), 
#              rep("adults", num_per_list*2), rep("elderly", num_per_list*2), 
#              rep("twentyplus", num_per_list*2))
#   temp <-  c(rep("constant", num_per_list), rep("var", num_per_list))
#   variable <- c(rep(temp, num_strategies))
#   
#   baseline_cases <- compute_total_cases(list_all$`0`)
#   baseline_cases_var <- compute_total_cases(list_all_var$`0`)
#   
#   temp <- c(rep(baseline_cases, num_per_list), rep(baseline_cases_var, num_per_list))
#   baseline_cases <- c(rep(temp, num_strategies))
#   
#   reduction <- (1-(total_cases/baseline_cases))*100
#   
#   df <- data.frame(vax_avail, strat, reduction, variable)
# }

get_reduction_in_cases_df_new = function(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  total_cases <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total_cases[count] <- compute_total_cases_new(i)
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
  
  baseline_cases <- compute_total_cases_new(list_all$`0`)
  baseline_cases_var <- compute_total_cases_new(list_all_var$`0`)
  
  temp <- c(rep(baseline_cases, num_per_list), rep(baseline_cases_var, num_per_list))
  baseline_cases <- c(rep(temp, num_strategies))
  
  reduction <- (1-(total_cases/baseline_cases))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

get_reduction_in_cases_df_single = function(this_list){
  total_cases <- c(NA)
  
  count <- 1
  for (i in this_list){
    total_cases[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  
  num_strategies <- 1
  vax_avail <- c(rep(seq(0, 50, by = 1)))
  num_per_list <- 51
  strat <- c(rep("optimal", num_per_list))
  variable <-  c(rep("constant", num_per_list))
  
  baseline_cases <- compute_total_cases_new(this_list$`0`)
  baseline_cases<- c(rep(baseline_cases, num_per_list))
  
  reduction <- (1-(total_cases/baseline_cases))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

get_reduction_in_cases_df_NTB = function(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  total_cases <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_cases[count] <- compute_total_cases_NTB(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total_cases[count] <- compute_total_cases_NTB(i)
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
  
  baseline_cases <- compute_total_cases_NTB(list_all$`0`)
  baseline_cases_var <- compute_total_cases_NTB(list_all_var$`0`)
  
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

get_reduction_in_deaths_df_new = function(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  total_deaths <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total_deaths[count] <- compute_total_deaths_new(i)
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
  
  baseline_deaths <- compute_total_deaths_new(list_all$`0`)
  baseline_deaths_var <- compute_total_deaths_new(list_all_var$`0`)
  
  temp <- c(rep(baseline_deaths, num_per_list), rep(baseline_deaths_var, num_per_list))
  baseline_deaths <- c(rep(temp, num_strategies))
  
  reduction <- (1 - (total_deaths/baseline_deaths))*100
  
  # df <- data.frame(vax_avail, strat, total_deaths, reduction, variable)
  df <- data.frame(vax_avail, strat, reduction, variable)
}

# get_reduction_in_deaths_df_novar = function(){
#   total_deaths <- rep(NA, 255)
#   
#   count <- 1
#   for (i in list_all){
#     total_deaths[count] <- compute_total_deaths(i)
#     count <- count + 1
#   }
#   for (i in list_kids){
#     total_deaths[count] <- compute_total_deaths(i)
#     count <- count + 1
#   }
#   for (i in list_adults){
#     total_deaths[count] <- compute_total_deaths(i)
#     count <- count + 1
#   }
#   for (i in list_elderly){
#     total_deaths[count] <- compute_total_deaths(i)
#     count <- count + 1
#   }
#   for (i in list_twentyplus){
#     total_deaths[count] <- compute_total_deaths(i)
#     count <- count + 1
#   }
#   
#   num_strategies <- 5
#   vax_avail <- c(rep(seq(0, 50, by = 1), num_strategies))
#   num_per_list <- 51
#   strat <- c(rep("all", num_per_list), rep("kids", num_per_list), 
#              rep("adults", num_per_list), rep("elderly", num_per_list), 
#              rep("twentyplus", num_per_list))
#   temp <-  c(rep("constant", num_per_list))
#   variable <- c(rep(temp, num_strategies))
#   
#   baseline_deaths <- compute_total_deaths(list_all$`0`)
#   
#   temp <- c(rep(baseline_deaths, num_per_list))
#   baseline_deaths <- c(rep(temp, num_strategies))
#   
#   reduction_in_deaths <- (1 - (total_deaths/baseline_deaths))*100
#   
#   df <- data.frame(vax_avail, strat, total_deaths, reduction_in_deaths, variable)
# }

get_reduction_in_deaths_df_single = function(this_list){
  total_deaths <- c(NA)
  
  count <- 1
  for (i in this_list){
    total_deaths[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  
  num_strategies <- 1
  vax_avail <- c(rep(seq(0, 50, by = 1)))
  num_per_list <- 51
  strat <- c(rep("optimal", num_per_list))
  variable <-  c(rep("constant", num_per_list))
  
  baseline_deaths <- compute_total_deaths_new(this_list$`0`)
  baseline_deaths <- c(rep(baseline_deaths, num_per_list))
  
  reduction <- (1-(total_deaths/baseline_deaths))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
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

get_reduction_in_YLL_df_new = function(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  total_YLL <- rep(NA, 510)
  
  count <- 1
  for (i in list_all){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_all_var){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_kids){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_kids_var){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_adults){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_adults_var){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_elderly){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_elderly_var){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  for (i in list_twentyplus_var){
    total_YLL[count] <- compute_total_YLL_new(i)
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
  
  baseline_YLL <- compute_total_YLL_new(list_all$`0`)
  baseline_YLL_var <- compute_total_YLL_new(list_all_var$`0`)
  
  temp <- c(rep(baseline_YLL, num_per_list), rep(baseline_YLL_var, num_per_list))
  baseline_YLL <- c(rep(temp, num_strategies))
  
  reduction <- (1 - (total_YLL/baseline_YLL))*100
  
  # df <- data.frame(vax_avail, strat, total_YLL, reduction, variable)
  df <- data.frame(vax_avail, strat, reduction, variable)
}

get_reduction_in_YLL_df_single = function(this_list){
  total_YLL <- c(NA)
  
  count <- 1
  for (i in this_list){
    total_YLL[count] <- compute_total_YLL_new(i)
    count <- count + 1
  }
  
  num_strategies <- 1
  vax_avail <- c(rep(seq(0, 50, by = 1)))
  num_per_list <- 51
  strat <- c(rep("optimal", num_per_list))
  variable <-  c(rep("constant", num_per_list))
  
  baseline_YLL <- compute_total_YLL_new(this_list$`0`)
  baseline_YLL <- c(rep(baseline_YLL, num_per_list))
  
  reduction <- (1-(total_YLL/baseline_YLL))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
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
                         Recovered = apply(df[29:37], 1, sum), 
                         Vaccinated = apply(df[38:46], 1, sum))
  
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
  theme_set(theme_minimal(base_size = 12))
  
  # string_to_add <- paste0("Total infected: ", round(compute_total_cases(df),3), "%")
  
  # Create a text
  # grob <- grobTree(textGrob(string_to_add, x=0.5,  y=0.95, hjust=0,
  #                           gp=gpar(col="white", fontsize=13, fontface="bold")))
  
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
      scale_y_continuous(expand = c(0,0), limit = c(0,15), breaks=c(0, 10, 20))
  } else if (compartment == "R") {
    p <- p + ylab("Percent Recovered") + ylim(0,100)
  } else if (compartment == "V") {
    p <- p + ylab("Percent Vaccinated") + ylim(0,100)
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

# Plot over vaccine avaliablity ----
plot_over_vax_avail = function(outcome, var_name, list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  
  theme_set(theme_minimal(base_size = 12))
  
  percent_kids <- ceiling(sum(age_demo[1:2])*100)
  percent_elderly <- ceiling(sum(age_demo[7:9])*100)
  percent_adults <- ceiling(sum(age_demo[3:5])*100)
  percent_twentyplus <- ceiling(sum(age_demo[3:9])*100)
  
  # get dataframe for specific outcome
  if (outcome == "cases"){
    df <- get_reduction_in_cases_df(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)
  } else if (outcome == "deaths"){
    df <- get_reduction_in_deaths_df(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)
    #df <- rbind(df[df$variable == "var", ], optimal_df_C[optimal_df_C$variable == "var", ])
  } else if (outcome == "YLL"){
    df <- get_reduction_in_YLL_df(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)
  }
  
  # lines stop when age groups are completely vaccinated
  # if (percent_kids < 50){df[df$strat == "kids" & df$vax_avail > percent_kids, ]$reduction = NA}
  # if (percent_elderly < 50) {df[df$strat == "elderly" & df$vax_avail > percent_elderly, ]$reduction = NA}
  # if (percent_adults < 50) {df[df$strat == "adults" & df$vax_avail > percent_adults, ]$reduction = NA}
  # if (percent_twentyplus < 50) {df[df$strat == "twentyplus" & df$vax_avail > percent_twentyplus, ]$reduction = NA}
  # 
  
  # reallocate to everyone after age groups are completely vaccinated according to the strat
  points_x <- c(percent_kids, percent_elderly, percent_adults)
  
  points_y <- c(df[df$strat == "kids" & df$vax_avail == floor(percent_kids) & df$variable == "constant", ]$reduction,
                df[df$strat == "elderly" & df$vax_avail == floor(percent_elderly) & df$variable == "constant", ]$reduction,
                df[df$strat == "adults" & df$vax_avail == floor(percent_adults) & df$variable == "constant", ]$reduction)
  points_y_var <- c(df[df$strat == "kids" & df$vax_avail == floor(percent_kids) & df$variable == "var", ]$reduction,
                    df[df$strat == "elderly" & df$vax_avail == floor(percent_elderly) & df$variable == "var", ]$reduction,
                    df[df$strat == "adults" & df$vax_avail == floor(percent_adults) & df$variable == "var", ]$reduction)
  points <- data.frame(points_x, points_y)
  
  # plot
  # for fig 2
  df_var <- df[df$variable == "var", ]
  df_const <- df[df$variable == "constant", ]
  p <- ggplot() +
    # geom_line(aes(x = df$vax_avail, y = df$reduction, col = df$strat, linetype = df$variable), size = 2, alpha = 0.9) +
    geom_line(aes(x = df_var$vax_avail, y = df_var$reduction, col = df_var$strat), linetype = "dashed", size = 2, alpha = 0.9) +
    geom_line(aes(x = df_const$vax_avail, y = df_const$reduction, col = df_const$strat), linetype = "solid", size = 1.2, alpha = 0.9) +
    xlab("Total vaccine supply (% of pop)") +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Adults 20-49", "All Ages", "Adults 60+", 
                                   "Under 20", "Adults 20+")) +
    scale_linetype_discrete(name = var_name,
                            labels = c("Constant", "Age-dependent")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100.2)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    #ggtitle(paste0(country)) +
    theme(legend.position = "none",
          #axis.title.x =element_blank(),
          #axis.text.x  = element_blank(),
          #axis.text.y  = element_blank(),
          #axis.title.y = element_blank()
    ) +
    theme(legend.text = element_text(size=22),
          legend.title = element_text(size =22),
          title = element_text(size = 22)) +
    guides(colour = guide_legend(override.aes = list(size=6)))
  
  p <- p + geom_point(aes(x = points_x, y = points_y), size = 2) + 
    geom_point(aes(x = points_x, y = points_y_var), size = 2) 
  
  if (outcome == "cases"){ p <- p + ylab("Reduction in infections (%)")}
  else if (outcome == "deaths"){ p <- p + ylab("Reduction in deaths (%)")}
  else if (outcome == "YLL"){ p <- p + ylab("Reduction in YLL (%)")}
  return(p)
}

ve_legend = function(){
  
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
  
  for (i in seq(0, 1, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim_new(C, j, "all", 1, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim_new(C, j, "kids", 1, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim_new(C, j, "adults", 1, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", 1, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", 1, v_e_type, this_v_e)
    list_all_var[[paste0(i)]] <- run_sim_new(C, j, "all", 1, v_e_type, v_e_var)
    list_kids_var[[paste0(i)]] <- run_sim_new(C, j, "kids", 1, v_e_type, v_e_var)
    list_adults_var[[paste0(i)]] <- run_sim_new(C, j, "adults", 1, v_e_type, v_e_var)
    list_elderly_var[[paste0(i)]] <- run_sim_new(C, j, "elderly", 1, v_e_type, v_e_var)
    list_twentyplus_var[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", 1, v_e_type, v_e_var)
  }

  df <- get_reduction_in_cases_df(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)

  p <- ggplot(df, aes(x = vax_avail, y = reduction, fill = strat)) +
    geom_line(size = 1, alpha = 1, aes(linetype = variable))+
    xlab("Total vaccine supply (% of pop)") +
    scale_linetype_discrete(name = "Efficacy",
                            labels = c("Constant", "Variable")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100.2)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    theme(legend.text = element_text(size=8),
          legend.title = element_text(size=8),
          legend.key.size = unit(0.85, "lines"),
          legend.direction = "vertical",
          legend.justification = c(1, 1),
          legend.margin=margin(t=1, r=0, b=-0.5, l=0, unit="cm")) +
    guides(shape = guide_legend(override.aes = list(size=0.5)))
  
  leg <- get_legend(p)
  plot(leg)
  return(leg)
}

color_legend = function(){
  
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")

  for (i in seq(0, 1, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim_new(C, j, "all", 1, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim_new(C, j, "kids", 1, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim_new(C, j, "adults", 1, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", 1, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", 1, v_e_type, this_v_e)
  }
  
  df <- get_reduction_in_cases_df_new(list_all , list_kids , list_adults , list_elderly , list_twentyplus )
  
  p <- ggplot(df, aes(x = vax_avail, y = reduction, col = strat)) +
    geom_line(size = 1, alpha = 1)+
    scale_color_brewer(palette = "Dark2", name = "Prioritization\nstrategy",
                       labels =  c("20-49", "All Ages", "60+", 
                                   "< 20", "20+")) +
    theme(legend.text = element_text(size=10),
          legend.title = element_text(size=10),
          legend.key.size = unit(0.1, "lines"),
          legend.justification = c(1, 1),
          legend.direction = "vertical",
          legend.margin=margin(t=0, r=0.3, b=0, l=0, unit="cm")) +
    guides(colour = guide_legend(override.aes = list(size=2)))
  
  leg <- get_legend(p)

  return(leg)
}

sero_legend = function(){
  
  list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")
  list_all_var <- list_kids_var <- list_adults_var <- list_elderly_var <- list_twentyplus_var <- vector(mode = "list")
  
  v_e_var <- get_v_e(0.5, 1, 50)
  for (i in seq(0, 1, by = 1)){
    j <- i/100
    list_all[[paste0(i)]] <- run_sim_new(C, j, "all", 1, v_e_type, this_v_e)
    list_kids[[paste0(i)]] <- run_sim_new(C, j, "kids", 1, v_e_type, this_v_e)
    list_adults[[paste0(i)]] <- run_sim_new(C, j, "adults", 1, v_e_type, this_v_e)
    list_elderly[[paste0(i)]] <- run_sim_new(C, j, "elderly", 1, v_e_type, this_v_e)
    list_twentyplus[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", 1, v_e_type, this_v_e)
    list_all_var[[paste0(i)]] <- run_sim_new(C, j, "all", 1, v_e_type, v_e_var)
    list_kids_var[[paste0(i)]] <- run_sim_new(C, j, "kids", 1, v_e_type, v_e_var)
    list_adults_var[[paste0(i)]] <- run_sim_new(C, j, "adults", 1, v_e_type, v_e_var)
    list_elderly_var[[paste0(i)]] <- run_sim_new(C, j, "elderly", 1, v_e_type, v_e_var)
    list_twentyplus_var[[paste0(i)]] <- run_sim_new(C, j, "twentyplus", 1, v_e_type, v_e_var)
  }
  
  df <- get_reduction_in_cases_df_new(list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var)
  
  p <- ggplot(df, aes(x = vax_avail, y = reduction, fill = strat)) +
    geom_line(size = 1, alpha = 1, aes(linetype = variable))+
    xlab("Total vaccine supply (% of pop)") +
    scale_linetype_discrete(name = "Distribute",
                            labels = c("w/o POC", "with POC")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100.2)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    theme(legend.text = element_text(size=10),
          legend.title = element_text(size =10),
          legend.key.size = unit(1, "lines"),
          legend.direction = "vertical",
          legend.justification = c(1, 1),
          legend.margin=margin(t=0, r=0.1, b=-0.5, l=0, unit="cm"))+
    guides(shape = guide_legend(override.aes = list(size=5)))
  
  leg <- get_legend(p)
  return(leg)
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

plot_over_vax_avail_new = function(outcome, var_name, list_all_var, list_kids_var, list_adults_var, list_elderly_var, list_twentyplus_var){
  library(ggplot2)
  theme_set(theme_minimal(base_size = 12))
  
  x_adults_switch     <- when_strat_switch(list_adults, 3:5)
  x_kids_switch       <- when_strat_switch(list_kids, 1:2)
  x_elderly_switch    <- when_strat_switch(list_elderly, 7:9)
  x_twentyplus_switch <- when_strat_switch(list_twentyplus, 3:9)
  
  x_adults_switch_var     <- when_strat_switch(list_adults_var, 3:5)
  x_kids_switch_var       <- when_strat_switch(list_kids_var, 1:2)
  x_elderly_switch_var    <- when_strat_switch(list_elderly_var, 7:9)
  x_twentyplus_switch_var <- when_strat_switch(list_twentyplus_var, 3:9)
  
  # get dataframe for specific outcome
  if (outcome == "cases"){
    df <- get_reduction_in_cases_df_new(list_all_var, list_kids_var, list_adults_var, list_elderly_var, 
                                        list_twentyplus_var)
    # df <- rbind(df, optimal_df_cases)
  } else if (outcome == "deaths"){
    df <- get_reduction_in_deaths_df_new(list_all_var, list_kids_var, list_adults_var, list_elderly_var, 
                                         list_twentyplus_var)
    # df <- rbind(df, optimal_df_deaths)
  } else if (outcome == "YLL"){
    df <- get_reduction_in_YLL_df_new(list_all_var, list_kids_var, list_adults_var, list_elderly_var, 
                                      list_twentyplus_var)
    # df <- rbind(df, optimal_df_YLL)
  }
  
  # reallocate to everyone after age groups are completely vaccinated according to the strat
  points_x <- c(x_kids_switch, x_elderly_switch, x_adults_switch)
  points_x_var <- c(x_kids_switch_var, x_elderly_switch_var, x_adults_switch_var)
  
  points_y <- c(df[df$strat == "kids" & df$vax_avail == x_kids_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "elderly" & df$vax_avail == x_elderly_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "adults" & df$vax_avail == x_adults_switch & df$variable == "constant", ]$reduction)
  points_y_var <- c(df[df$strat == "kids" & df$vax_avail == x_kids_switch_var & df$variable == "var", ]$reduction,
                    df[df$strat == "elderly" & df$vax_avail == x_elderly_switch_var & df$variable == "var", ]$reduction,
                    df[df$strat == "adults" & df$vax_avail == x_adults_switch_var & df$variable == "var", ]$reduction)
  
  if (x_twentyplus_switch > 0){
    points_x <- c(points_x, x_twentyplus_switch)
    points_y <- c(points_y, df[df$strat == "twentyplus" & df$vax_avail == x_twentyplus_switch & df$variable == "constant", ]$reduction)
    points_y_var <- c(points_y_var, df[df$strat == "twentyplus" & df$vax_avail == x_twentyplus_switch & df$variable == "var", ]$reduction)
  }
  
  # plot
  # for fig 2
  df_var <- df[df$variable == "var", ]
  df_const <- df[df$variable == "constant", ]
  p <- ggplot() +
    # geom_line(aes(x = df_var$vax_avail, y = df_var$reduction, col = df_var$strat), 
    #          linetype = "dashed", size = 1, alpha = 0.9) +
    geom_line(aes(x = df_const$vax_avail, y = df_const$reduction, col = df_const$strat), 
              linetype = "solid", size = 1, alpha = 0.9) +
    geom_line(aes(x = df_var[df_var$strat == "optimal",]$vax_avail, 
                  y = df_var[df_var$strat == "optimal",]$reduction, 
                  col = df_var[df_var$strat == "optimal",]$strat),
              linetype = "solid", size = 1, alpha = 0.9) +
    xlab("Total vaccine supply (% of population)") +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Adults 20-49", "All Ages", "Adults 60+", 
                                   "Under 20", "Adults 20+")) +
    scale_linetype_discrete(name = df$var_name,
                            labels = c("Constant", "Age-dependent")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    coord_fixed(0.5*4/5) +
    theme(legend.position = "none") +
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  p <- p + geom_point(aes(x = points_x, y = points_y), size = 1.5) + 
    geom_point(aes(x = points_x_var, y = points_y_var), size = 1.5) 
  
  if (outcome == "cases"){ p <- p + ylab("Reduction in\ninfections (%)")}
  else if (outcome == "deaths"){ p <- p + ylab("Reduction in\ndeaths (%)")}
  else if (outcome == "YLL"){ p <- p + ylab("Reduction in\nYLL (%)")}
  return(p)
}

plot_over_vax_avail_NTB = function(outcome, var_name, list_all_var, list_kids_var, 
                                   list_adults_var, list_elderly_var, list_twentyplus_var){
  library(ggplot2)
  theme_set(theme_minimal(base_size = 12))
  
  x_adults_switch     <- when_strat_switch(list_adults, 3:5)
  x_kids_switch       <- when_strat_switch(list_kids, 1:2)
  x_elderly_switch    <- when_strat_switch(list_elderly, 7:9)
  x_twentyplus_switch <- when_strat_switch(list_twentyplus, 3:9)
  
  x_adults_switch_var     <- when_strat_switch(list_adults_var, 3:5)
  x_kids_switch_var       <- when_strat_switch(list_kids_var, 1:2)
  x_elderly_switch_var    <- when_strat_switch(list_elderly_var, 7:9)
  x_twentyplus_switch_var <- when_strat_switch(list_twentyplus_var, 3:9)
  
  # get dataframe for specific outcome
  if (outcome == "cases"){
    df <- get_reduction_in_cases_df_NTB(list_all_var, list_kids_var, list_adults_var, list_elderly_var, 
                                        list_twentyplus_var)
  } else if (outcome == "deaths"){
    df <- get_reduction_in_deaths_df_NTB(list_all_var, list_kids_var, list_adults_var, list_elderly_var, 
                                         list_twentyplus_var)
    #df <- rbind(df[df$variable == "var", ], optimal_df_C[optimal_df_C$variable == "var", ])
  } else if (outcome == "YLL"){
    df <- get_reduction_in_YLL_df_NTB(list_all_var, list_kids_var, list_adults_var, list_elderly_var, 
                                      list_twentyplus_var)
  }
  
  # reallocate to everyone after age groups are completely vaccinated according to the strat
  points_x <- c(x_kids_switch, x_elderly_switch, x_adults_switch)
  points_x_var <- c(x_kids_switch_var, x_elderly_switch_var, x_adults_switch_var)
  
  points_y <- c(df[df$strat == "kids" & df$vax_avail == x_kids_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "elderly" & df$vax_avail == x_elderly_switch & df$variable == "constant", ]$reduction,
                df[df$strat == "adults" & df$vax_avail == x_adults_switch & df$variable == "constant", ]$reduction)
  points_y_var <- c(df[df$strat == "kids" & df$vax_avail == x_kids_switch_var & df$variable == "var", ]$reduction,
                    df[df$strat == "elderly" & df$vax_avail == x_elderly_switch_var & df$variable == "var", ]$reduction,
                    df[df$strat == "adults" & df$vax_avail == x_adults_switch_var & df$variable == "var", ]$reduction)
  
  if (x_twentyplus_switch > 0){
    points_x <- c(points_x, x_twentyplus_switch)
    points_y <- c(points_y, df[df$strat == "twentyplus" & df$vax_avail == x_twentyplus_switch & df$variable == "constant", ]$reduction)
    points_y_var <- c(points_y_var, df[df$strat == "twentyplus" & df$vax_avail == x_twentyplus_switch & df$variable == "var", ]$reduction)
  }
  
  # plot
  # for fig 2
  df_var <- df[df$variable == "var", ]
  df_const <- df[df$variable == "constant", ]
  p <- ggplot() +
    geom_line(aes(x = df_var$vax_avail, y = df_var$reduction, col = df_var$strat), 
              linetype = "dashed", size = 1, alpha = 0.9) +
    geom_line(aes(x = df_const$vax_avail, y = df_const$reduction, col = df_const$strat), 
              linetype = "solid", size = 0.5, alpha = 0.9) +
    xlab("Total vaccine supply (% of population)") +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Adults 20-49", "All Ages", "Adults 60+", 
                                   "Under 20", "Adults 20+")) +
    scale_linetype_discrete(name = df$var_name,
                            labels = c("Constant", "Age-dependent")) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50)) +#, breaks = c(0,25,50)) +
    coord_fixed(0.5*4/5) +
    theme(legend.position = "none") +
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  p <- p + geom_point(aes(x = points_x, y = points_y), size = 1.5) + 
    geom_point(aes(x = points_x_var, y = points_y_var), size = 1.5) 
  
  if (outcome == "cases"){ p <- p + ylab("Reduction in\ninfections (%)")}
  else if (outcome == "deaths"){ p <- p + ylab("Reduction in\ndeaths (%)")}
  else if (outcome == "YLL"){ p <- p + ylab("Reduction in\nYLL (%)")}
  return(p)
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
                   axis.text.x = element_text(size = 17), 
                   axis.text.y = element_text(size = 17),
                   plot.title = element_blank(),
                   legend.position = "none") 
  }
  
  for (i in 1:100){
    p <- p + geom_line(data = df[df$prev_num == paste0(i),],
                       aes(x = vax_avail, y = reduction, col = strat),
                       size = 1.1, alpha = 0.25)
  }
  
  return(p)
} 

plot_serotesting_over_vax_avail = function(outcome){
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
    
    reduction <- (1-(total_cases/baseline_cases))*100
  }
  else if (outcome == "deaths"){
    baseline_deaths_notest <- compute_total_deaths(list_all_sero_notest$`0`)
    baseline_deaths_test <- compute_total_deaths(list_all_sero_test$`0`)
    
    temp <- c(rep(baseline_deaths_notest, num_per_list), rep(baseline_deaths_test, num_per_list))
    baseline_deaths <- c(rep(temp, num_strategies))
    
    reduction <- (1-(total_deaths/baseline_deaths))*100
  }
  df <- data.frame(vax_avail, strat, reduction, variable)
  # for Belgium sero
  df[df$strat == "kids" & vax_avail > 23 & variable == "notest", ]$reduction = NA
  df[df$strat == "kids" & vax_avail > 22 & variable == "test", ]$reduction = NA
  df[df$strat == "elderly" & vax_avail > 26 & variable == "notest", ]$reduction = NA
  df[df$strat == "elderly" & vax_avail > 25 & variable == "test", ]$reduction = NA
  df[df$strat == "adults" & vax_avail > 38 & variable == "notest", ]$reduction = NA
  df[df$strat == "adults" & vax_avail > 36 & variable == "test", ]$reduction = NA
  
  # ## for NY sero
  # df[df$strat == "kids" & vax_avail > percent_kids & variable == "notest", ]$reduction_in_cases = NA
  # df[df$strat == "kids" & vax_avail > 18 & variable == "test", ]$reduction_in_cases = NA
  # df[df$strat == "elderly" & vax_avail > percent_elderly & variable == "notest", ]$reduction_in_cases = NA
  # df[df$strat == "elderly" & vax_avail > 18 & variable == "test", ]$reduction_in_cases = NA
  # df[df$strat == "adults" & vax_avail > percent_adults & variable == "notest", ]$reduction_in_cases = NA
  # df[df$strat == "adults" & vax_avail > 30 & variable == "test", ]$reduction_in_cases = NA
  # df[df$strat == "all" & vax_avail > 37 & variable == "test", ]$reduction_in_cases = NA
  
  ggplot(df, aes(x = vax_avail, y = reduction, col = strat)) +
    geom_line(aes(linetype = variable), size = 1.2, alpha = 0.9) +
    xlab("Total vaccine supply (% of pop)") +
    ylab("Reduction in infections (%)") +
    scale_linetype_discrete(name = "Scenario",
                            labels = c("w/o sero tests", "w/ sero tests"))+
    theme(
      #legend.text = element_text(size=22),
      #legend.title = element_text(size=22),
      legend.position = "none",
      #axis.title.x = element_blank(),
      axis.text.y  = element_blank(),
      axis.text.x = element_text(size = 17),
      #axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank()) +
    scale_y_continuous(expand = c(0,0), limit = c(0, 100)) +
    scale_x_continuous(expand = c(0,0), limit = c(0, 50), breaks = c(0,10, 20, 30, 40, 50)) +
    scale_color_brewer(palette = "Dark2", name = "Allocation Strategy",
                       labels =  c("Adults 20-49", "All Ages", "Adults 60+", 
                                   "Under 20", "Adults 20+", "Optimal")) +
    guides(colour = guide_legend(override.aes = list(size=6)))
}

# Other plotting fn ----
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
                                labels = seq(1.3, 2.6, by = 0.2))
  } else if (param == "ve"){
    p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 2),
                                labels = c(" 30", " 50", " 70", " 90"))
  } else if (param == "country") {
    p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 1),
                              labels = c("BEL", "USA", "IND", "ESP", "ZWE", "BRA", "CHN", "ZAF", "POL"))
  } else if (param == "rollout") {
    p <- p + scale_y_discrete(expand = c(0,0), breaks = seq(1, num, by = 1),
                              labels = seq(0.25, 2, by = 0.25))
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

# OLD PLOTTING CODE FN - needs to be updated
# barplot_at_finalT = function(solution_df, IFR, strategy, y, nvax = 0){
#   if (strategy == "no vax"){
#     vaccinated <- rep(0,num_groups) * 100
#     this_color <- "#808080"
#     plot_title <- "No Vaccines"
#   } 
#   else {
#     if (strategy == "all"){
#       groups <- 1:9
#       this_color <- col_all
#       plot_title <- "All Ages"
#       sub_title <- ""
#     } 
#     else if (strategy == "kids"){
#       groups <- 1:2
#       this_color <- col_kids
#       plot_title <- "Children (0-19)"
#       sub_title <- ""
#     } 
#     else if (strategy == "adults"){
#       groups <- 3:5
#       this_color <- col_youngadults
#       plot_title <- "Young Adults (20-49)"
#       sub_title <- "(20-49)"
#     } 
#     else if (strategy == "elderly"){
#       groups <- 7:9
#       this_color <- col_elderly
#       plot_title <- "Elderly (60+)"
#       sub_title <- "(60+)"
#     }
#     else if (strategy == "20+"){
#       groups <- 3:9
#       this_color <- col_adults
#       plot_title <- "Adults (20+)"
#       sub_title <- "(20+)"
#     }
#     people_to_vax <- sum(N_i[groups])
#     
#     vax_proportion <- rep(0, num_groups)
#     vax_proportion[groups] <- N_i[groups]/people_to_vax
#     
#     vax_distribution <- nvax*vax_proportion
#     vax_distribution[vax_distribution > N_i] <- N_i[vax_distribution > N_i]
#     vaccinated <- vax_distribution/N_i * 100
#   }
#   
#   infections <- rep(0,num_groups)
#   deaths <- rep(0, num_groups)
#   final_S <- rep(0, num_groups)
#   
#   R_index <- 29 # col number for R1
#   
#   for (i in 1:num_groups) {
#     deaths[i] <- max(solution_df[R_index])*IFR[i]/N_i[i]*100
#     infections[i] <- max(solution_df[R_index]) # total infected (including dead)
#     # infections[i] <- max(solution_df[R_index])*(1-IFR[i])/N_i[i] # total infected (excluding dead)
#     final_S[i] <- solution_df[(max(t) + 1), (i + 1)]/N_i[i]
#     R_index <- R_index + 1
#   }
#   
#   #age_groups <- c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
#   age_groups <- c("0","10", "20", "30", "40", "50", "60", "70", "80")
#   
#   df <- data.frame(age_groups, final_S, deaths, infections, vaccinated)
#   
#   new_df <- df %>%
#     # select(age_groups, infections, deaths, vaccinated, final_S) %>%
#     select(age_groups, all_of(y)) %>%
#     gather(key = "type", value = "num", -age_groups)
#   
#   theme_set(theme_classic(base_size = 20))
#   
#   # plot
#   p <- ggplot(new_df, aes(fill=type, y=num, x=age_groups)) + 
#     geom_bar(position="stack", stat="identity", fill = this_color) + 
#     xlab("Age group") +
#     scale_x_discrete(breaks=c("0","10", "20", "30", "40", "50", "60","70", "80"), 
#                      labels=c("0","", "20", "", "40", "", "60","", "80")) + 
#     scale_y_continuous(expand = c(0,0), limit = c(0, 100), breaks = c(0, 50, 100)) +
#     theme(plot.subtitle = element_text(size = 9, face = "italic", hjust = 1),
#           plot.title = element_text(size = 17),
#           axis.title.x =element_blank(), 
#           axis.title.y = element_blank())
#   
#   if (y == "infections") {
#     string_to_add <- paste0("Total infected: ", round(compute_total_cases(solution_df),2), "%")
#     
#     p <- p + 
#       #ylim(0, 100) + 
#       ylab("Infected (%)") + 
#       # labs(title = plot_title, subtitle = string_to_add)
#       labs(subtitle = string_to_add)
#   } 
#   else if (y == "deaths") {
#     string_to_add <- paste0("Total deaths: ", round(compute_total_deaths(solution_df),2), "%")
#     
#     p <- p + 
#       ylim(0, 8) + 
#       ylab("Deaths (%)") + 
#       # labs(title = plot_title, subtitle = string_to_add)
#       labs(subtitle = string_to_add)
#   } 
#   else if (y == "vaccinated"){
#     p <- p + 
#       #ylim(0, 100) + 
#       ylab("Vaccinated (%)") 
#     #labs(title = plot_title)
#   }
#   p # output ggplot
# }
# plot_totalinfections_sero =function(df_baseline, df_sero_notest, df_sero_test, compartment){
#   # INPUTS:
#   # df: simulation to plot
#   # compartment: character of compartment to plot i.e. "I", "R"
#   #
#   # OUTPUT:
#   # p: ggplot object of plot
#   
#   # make new df for plotting
#   if (compartment == "I"){
#     total <- data.frame(time = 0:150,
#                         baseline = apply(df_baseline[20:28], 1, sum), 
#                         sero_notest = apply(df_sero_notest[20:28], 1, sum),
#                         sero_test = apply(df_sero_test[20:28], 1, sum))
#   } else if (compartment == "R"){
#     total <- data.frame(time = 0:150,
#                         baseline = apply(df_baseline[29:37], 1, sum), 
#                         sero_notest = apply(df_sero_notest[29:37], 1, sum),
#                         sero_test = apply(df_sero_test[29:37], 1, sum))
#   }
#   
#   new_df <- total %>%
#     select(time, baseline, sero_notest, sero_test) %>%
#     gather(key = "scenario", value = "num", -time)
#   
#   # Add col for percent infected & recovered
#   new_df$percent <- new_df$num/pop_total*100
#   
#   # Plot
#   theme_set(theme_minimal(base_size = 20))
#   
#   # string_to_add <- paste0("Total infected: ", compute_total_cases(df), "%")
#   
#   # Create a text
#   # grob <- grobTree(textGrob(string_to_add, x=0.5,  y=0.95, hjust=0,
#   #                           gp=gpar(col="white", fontsize=13, fontface="bold")))
#   
#   # Plot
#   p <- ggplot(new_df, aes(x = time, y = percent)) + 
#     geom_line(aes(color = scenario), size = 2.5) +
#     xlab("Time") +
#     guides(colour = guide_legend(override.aes = list(size=6))) +
#   #   annotation_custom(grob) +
#     scale_color_manual(values = wes_palette("GrandBudapest1", n = 3), name = "Scenario", 
#                         labels =  c("Baseline\n(no seropositive)\n", "Seropositive\nno testing\n", "Seropositive\nwith testing"))
#   # 
#   if (compartment == "I") {
#     p <- p + ylab("Percent Infected") + ylim(0, 16)
#   } else if (compartment == "R") {
#     p <- p + ylab("Percent Recovered") + ylim(0,1)
#   }
# }