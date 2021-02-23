# Vaccine_Allocation_Project

Code layout: 

main.R - script to run the model for evaluating the 5 different strategies and most of the plot figures, enter the country of interest at the top of the script using the three letter country code. Organized into blocks by figure. Run each block at a time. 

main_heatmap.R, run_heatmap_R0.R, run_heatmap_rollout.R - script to run the model and produce corresponding heatmap plots to visualize the best strategy under specified parameter values

helper_functions.R - most of the functions for main.R, model related functions at the top, plotting functions at the bottom
                   
run_sim.R - model function to run the simulation run_sim (for all-or-nothing efficacy) or run_sim_NTB (for partially transmission blocking efficacy)

functions starting with "get_" - convert original data into .RData to be used in either main or main_optimal

set_up.R - loads packages, plotting functions and initializes all the default parameter settings

See https://vaxfirst.colorado.edu/ for a interactive calculator of this model and to investigate how different parameter values change the reduction in infections, deaths and years of life lost under each vaccine prioritization strategy. 



# License

Copyright (C) 2021, Kate Bubar

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (see LICENSE.md).  If not, see <https://www.gnu.org/licenses/>.
