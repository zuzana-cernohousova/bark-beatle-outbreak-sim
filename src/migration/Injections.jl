include("MigrationSimFunctions.jl")
using .MigrationSimFunctions

include("SpatialSimulations.jl")
using .SpatialSimulations

include("plotting.jl")
using .SpatialPlotting

include("Scenarios.jl")
using .Scenarios

using Plots
using Colors
using Random
using NamedTupleTools

Random.seed!(4269)

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
    gamma = 300.0,

    # from approximate values table
    lambda = 0.001, # only this value in table
    r = 0.1,        # only this value in table
    e = 10,         # lower end --> test in grid

    # aggregation function
    beta_function = aggreg_beta_mig,

    # my parameter
    a = 5
)

dimensions = (x = 10, y = 10)
orig_params_matrix = fill(orig_params, dimensions.x, dimensions.y)

S_init = Float64.(fill(1.15, dimensions.x, dimensions.y))
I_init = Float64.(fill(1.0, dimensions.x, dimensions.y))

# ----- start with increased I in position
for i in 1.1:0.1:1.4  
    S_init = Float64.(fill(i, dimensions.x, dimensions.y))
  
    playout_spatial_injection((S = S_init, I = I_init),
        orig_params_matrix,
        (x = 5, y = 4),
        (x = 5, y = 5),
        (S = 1.15, I = 2);
        show_gif = true, stop_t = 15,
        dir_path = "../../out/fig_8/")
end

# viz_for_init_at_pos((S = S_init, I = I_init),
# orig_params_matrix,
# [(x = 5, y = 5), (x =5, y = 4)],
# " without injection";
# show_gif = false, stop_t = 15,
# dir_path = "../../out/fig_7/")
