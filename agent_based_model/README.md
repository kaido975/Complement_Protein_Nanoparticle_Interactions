Open and Execute dependencies.jl (!!! For version control !!!)

To run the ABM
Open the terminal and execute
julia -p 4 run_distributed.jl

Data will be saved in #Number of agents..agents_data.csv (containing 500 simulations)

Use the R script to analyze the occupied sites (To reproduce density plots in Fig. 7e).

Run the simulation for different numbers of binding sites by changing line number 13 in run_distributed.jl
