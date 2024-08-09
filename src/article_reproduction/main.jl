include("SimulationFunctions.jl")
include("Simulations.jl")

using Plots
import PyPlot
using .SimulationFunctions
using .Simulations

orig_params = (
    # from fig 1
    g = 0.001,
    d = 0.003,
    m = 0.05,
    beta = 0.01, 
    mu = 2000, 
    K = 100,
    n = 10,

    # intermediate_resistance
    gamma = 300,

    # from approximate values table
    lambda = 0.001, # only this value in table
    r = 0.1,        # only this value in table
    e = 10,         # lower end

    # aggregation function
    beta_function = aggreg_beta,
)

# # simple simulation
# init = (
#     S = 1,
#     I = 1
# )

# Ss, Is = simple_simulation(orig_params, init, (dS = dS, dI = dI); delta_t = 0.1, stop_t = 5)
# display(plot([Ss, Is], xlabel = "Time", label=["S" "I"], linewidth=3)) # create plot of S and I

# run exploratory_analysis to visualize figure 3
ea_results = exploratory_analysis((S, I, nondim_params, orig_params) -> [dS(S, I, nondim_params), dI(S, I, nondim_params)], orig_params,
    "../../out/fig_3/",
    [
        [:beta_function, [aggreg_beta, uniform_beta]],
        [:gamma, [30, 300, 450]]
    ], 
    save_to_files = true
)

