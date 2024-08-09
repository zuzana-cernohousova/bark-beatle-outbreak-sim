module SpatialPlotting

export create_gif
export create_gamma_gradient
export create_line_graph_in_position
export visualize_simulation_results
export get_frame_in_step

using Plots
using Colors

function scale_matrices(matrices, new_min, new_max)
    # Find the global minimum and maximum
    global_min = minimum([minimum(matrix) for matrix in matrices])
    global_max = maximum([maximum(matrix) for matrix in matrices])
    
    # Define a nested function to scale a single matrix
    function scale_matrix(matrix, sub_min, sub_max)
        # Scale to [0, 1]
        scaled_matrix = (matrix .- global_min) ./ (global_max .- global_min)

        # invert if necessary
        if sub_min > sub_max
            sub_min, sub_max = sub_max, sub_min
            scaled_matrix = 1 .- scaled_matrix  # Reverse the scaling
        end

        # Scale to [new_min, new_max]
        return scaled_matrix .* (sub_max - sub_min) .+ sub_min
    end

    return [scale_matrix(matrix, new_min, new_max) for matrix in matrices]
end

function matrix_to_color(matrix1, matrix2)
    rows, cols = size(matrix1)
    img = Array{RGB{Float64}}(undef, rows, cols)
    
    for i in 1:rows
        for j in 1:cols
            hue = matrix2[i, j]
            intensity = matrix1[i, j]
            img[i, j] = HSV(hue, 1.0, intensity)
        end
    end
    
    return img
end

function create_gif(Ss, Is)
    a = Animation()
    
    max_S = maximum([maximum(matrix) for matrix in Ss])
    max_I = maximum([maximum(matrix) for matrix in Is])

    for i in 1:length(Ss)
        hm_S = heatmap(Ss[i], color=:linear_kbc_5_95_c73_n256, clim = (0, max_S), title= "Susceptible", aspect_ratio = 1.0)
        hm_I = heatmap(Is[i], color=:linear_kry_0_97_c73_n256, clim = (0, max_I), title = "Infested", aspect_ratio = 1.0)

        plt = plot(hm_I, hm_S, layout= (1,2), size = (800, 400))
        frame(a, plt)
    end
    return a
end

function get_frame_in_step(Ss, Is, step)
    max_S = maximum([maximum(matrix) for matrix in Ss])
    max_I = maximum([maximum(matrix) for matrix in Is])

    hm_S = heatmap(Ss[step], color=:linear_kbc_5_95_c73_n256, clim = (0, max_S), title= "Susceptible", aspect_ratio = 1.0)
    hm_I = heatmap(Is[step], color=:linear_kry_0_97_c73_n256, clim = (0, max_I), title = "Infested", aspect_ratio = 1.0)

    return plot(hm_I, hm_S, layout= (1,2), size = (800, 400))
end

function custom_gradient(x, minimum, maximum)
    minimum + (maximum-minimum)*(((-x*x) + 10*x)/25) 
end

function create_gamma_gradient(params_matrix, gamma_min, gamma_max)
    for x in 1:10
        for y in 1:10
            params_matrix[x,y] = merge(params_matrix[x,y], [:gamma => sqrt((custom_gradient(x, gamma_min, gamma_max)^2 + custom_gradient(y, gamma_min, gamma_max)^2)/2)])
        end
    end
end

function get_vector_in_position(matrices, pos)
    return [matrix[pos.x, pos.y] for matrix in matrices]
end

function create_line_graph_in_position(S_matrix, I_matrix, pos; title = "")
    plot([get_vector_in_position(S_matrix, pos), get_vector_in_position(I_matrix, pos)], xlabel = "Time", label=["S" "I"], linewidth=3, title=title, size = (600, 450), titlefontsize = 15) 
end

function visualize_simulation_results(results, observation_positions; show_gif = false, dir_path = "", comment = "")
    for op in observation_positions
        title = "position ["*string(op.x)*", "*string(op.y)*"]"* comment

        graph = create_line_graph_in_position(results..., op; title = title)

        if dir_path != ""
            savefig(dir_path * replace(string(title), ' ' => '_') *".png")
        else
            display(graph)
        end
    end

    if show_gif 
        if dir_path != ""
            gif(create_gif(results...), dir_path * replace(string(comment), ' ' => '_') *".gif", fps = 10)
        else
            gif(create_gif(results...), fps = 10)
        end
    end
end
end