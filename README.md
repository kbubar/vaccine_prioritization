# Vaccine_Allocation_Project

Code layout: 
main.R - script to run the model for evaluating the 5 different strategies and plot figures
main_optimal.R - script to run the optimization code for the model and plot some of the related optimal-related figures
helper_functions.R - most of the functions for main.R and main_optimal.R
                   - model related functions at the top, plotting functions at the bottom
runsim.R - model functions
         - includes running the basic simulation, nontransmission blocking simulation, WHOstrat (input who to vaccinate by percent instead of by strategy), and simulation to                  distribute vax after t = 0
