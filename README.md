# Vaccine_Allocation_Project

Code layout: 

main.R - script to run the model for evaluating the 5 different strategies and plot figures (includes a lot of code for old figures), enter the country of interest using the three letter country code

main_optimal.R - script to run the optimization code for the model and plot some of the related optimal-related figures

helper_functions.R - most of the functions for main.R and main_optimal.R, model related functions at the top, plotting functions at the bottom
                   
run_sim.R - model function
         - includes running the basic simulation, nontransmission blocking simulation, WHOstrat (input who to vaccinate by percent instead of by strategy), and simulation to                  distribute vax after t = 0

optimize_sim.R - function for computing optimal allocation

plot_optimal_over_vax.R - plotting script for computing the best allocation across multiple optimization simulations and plotting

functions starting with "get_" - convert original data into .RData to be used in either main or main_optimal
