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

    # my parameter for migration
    a = 10
)

dimensions = (x = 10, y = 10)
orig_params_matrix = fill(orig_params, dimensions.x, dimensions.y)

S_init = Float64.(fill(1.15, dimensions.x, dimensions.y))
I_init = Float64.(fill(1.0, dimensions.x, dimensions.y))

# ------ randomized gamma parameter with different a
randomize_parameter(orig_params_matrix, :gamma, 150.0, 450.0)

gammas = [t.gamma for t in orig_params_matrix]
display(heatmap(gammas, color=:linear_kgy_5_95_c69_n256, clim = (150, 450), title= "Gammas init", aspect_ratio = 1.0, size=(400,400)))
savefig("../../out/fig_10/gammas_init.png")

comment = "random gammas from 150 to 450, with a = " * string(orig_params.a)
random_gammas_result = spatial_simulation((S = S_init, I = I_init), orig_params_matrix; stop_t = 40)
visualize_simulation_results(random_gammas_result, [];dir_path = "../../out/fig_11/gifs/", show_gif = true, comment = comment)
savefig("../../out/fig_11/"*replace(comment, " " => "_")*".png")


# ------ one gamma set lower

# gammas = [30, 50, 100, 200]
# for g in gammas
#     playout_parameter_injection(
#         (S = S_init, I = I_init),
#         orig_params_matrix,
#         :gamma, g, (x = 5, y = 5),
#         [(x = 5, y = 5), (x = 5, y = 4), (x = 4, y = 4)];
#         stop_t = 15, show_gif = false,
#         dir_path = "../../out/fig_9/")
# end