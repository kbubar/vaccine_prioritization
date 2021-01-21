# _____________________________________________________________________
# IMPORT ----
# _____________________________________________________________________
library(tidyverse) 
library(deSolve) # ode solver
library(gridExtra)
library(RColorBrewer)
library(wesanderson)
library(egg)
library(foreach)
library(doParallel)
library(ggplot2)
library(gplots)

# ColorBrewer Dark2
col_youngadults = "#1B9E77" # teal green (20-49)
col_all = "#D95F02" # orange
col_elderly = "#7570B3" # purple (60+)
col_kids = "#E7298A" # magenta
col_adults = "#66A61E" # light green (20+ strategy)

tippingpointPal2 <- c("#E6E6E6", "#FFEDA0", "#FEB24C", "#F03B20", "#000000")
names(tippingpointPal2) <- c("None", "0-25%", "25-50%",
                             "50-100%", "NA")

colFill2 <- scale_fill_manual(name = "Tipping point", values = tippingpointPal2)

nolabels_theme <- theme(axis.title.x =element_blank(),
                        axis.text.x = element_blank(),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        plot.title = element_text(size = 12, face = "plain"),
                        legend.position = "none")
onlyx_theme <- theme(axis.title.y = element_blank(),
                     axis.text.y = element_blank(),
                     plot.title = element_text(size = 12, face = "plain"),
                     legend.position = "none")
onlyy_theme <- theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank(),
                     plot.title = element_text(size = 12, face = "plain"),
                     legend.position = "none")
alllabels_theme <- theme(plot.title = element_text(size = 12, face = "plain"),
                         legend.position = "none")

theme_set(theme_minimal(base_size = 12))

####
d_E <- 1/3 # incubation period (E -> I), ref: Davies
d_I <- 1/5 # recovery period (I -> R), ref: Davies

C          <- readRDS(paste0("C_", country, "_bytens_overall.RData"))

age_demo   <- readRDS(paste0("age_demographics_", country,".RData"))
pop_total  <- age_demo[10]
age_demo   <- age_demo[1:9]

N_i        <- pop_total*age_demo      
num_groups <- length(age_demo) # num age groups

#IFR        <- c(9.530595e-04, 3.196070e-03, 1.071797e-02, 3.594256e-02, 1.205328e-01, 
#                4.042049e-01, 1.355495e+00, 4.545632e+00, 1.524371e+01) # old IFR, Ref: Levin2020
IFR        <- c(9.807703e-04, 3.277686e-03, 1.095386e-02, 3.660727e-02, 1.223397e-01, 4.088531e-01, 1.366367e+00,
                4.566330e+00, 1.526045e+01) # updated IFR, Ref: Levin2020
IFR        <- IFR/100 # as decimal
YLL_vec    <- readRDS(paste0("yll_vec_", country, ".RData"))

u_var      <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74) # Ref: Davies2020

this_v_e   <- get_v_e(p = 0.9, y0 = 0.9, hinge_age = 50)
v_e_var    <- get_v_e(p = 0.5, y0 = 0.9, hinge_age = 50)
v_e_type   <- "aorn"

sero_none <- rep(0, 9) # no prior immunity
sero_CT <- c(0.039, 0.0382, 0.031, 0.031, 0.031, 0.037, 0.032, 0.027, 0.027)
sero_NY <- c(0.32, 0.3129, 0.249, 0.249, 0.264, 0.279, 0.2575, 0.2215, 0.207) # ref: NYC

this_sp <- 0.99
this_se <- 0.70

num_perday <- 0.002
list_all <- list_kids <- list_adults <- list_elderly <- list_twentyplus <- vector(mode = "list")

scale_115 <- scale_u_for_R0(u_var, C, 1.15)
scale_15 <- scale_u_for_R0(u_var, C, 1.5)
scale_26 <- scale_u_for_R0(u_var, C, 2.6)

C_115 <- C/scale_115
C_15  <- C/scale_15
C_26  <- C/scale_26