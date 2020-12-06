# Vaccine_Allocation_Project

Code layout: 

main.R - script to run the model for evaluating the 5 different strategies and plot figures (includes a lot of code for old figures), enter the country of interest using the three letter country code

helper_functions.R - most of the functions for main.R, model related functions at the top, plotting functions at the bottom
                   
run_sim.R - model function
         - includes running the basic simulation, nontransmission blocking simulation, WHOstrat (input who to vaccinate by percent instead of by strategy), and simulation to                  distribute vax after t = 0
         - the revised version of the model is called with run_sim_new

functions starting with "get_" - convert original data into .RData to be used in either main or main_optimal
