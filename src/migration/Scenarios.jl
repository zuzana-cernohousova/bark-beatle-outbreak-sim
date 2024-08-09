module Scenarios

using ..SpatialSimulations
using ..SpatialPlotting

export playout_spatial_injection
export playout_parameter_injection
export randomize_parameter
export viz_for_init_at_pos

function playout_spatial_injection(init, orig_params_matrix, observation_pos, injection_pos, inj_tuple; show_gif = false, stop_t = 10, dir_path = "")    
    orig_s = init.S[1,1]
    init.I[injection_pos...] = inj_tuple.I
    init.S[injection_pos...] = inj_tuple.S

    injection_result = spatial_simulation(init, orig_params_matrix, stop_t = stop_t)
    visualize_simulation_results(injection_result, (injection_pos, observation_pos); show_gif = show_gif, dir_path = dir_path, comment = " with init S = " * string(orig_s))
end

function viz_for_init_at_pos(init, orig_params_matrix, obs_pos, comment; show_gif = false, stop_t = 10, dir_path = "")
    res = spatial_simulation(init, orig_params_matrix, stop_t = stop_t)
    visualize_simulation_results(res, obs_pos; show_gif = show_gif, dir_path = dir_path, comment = comment)
end

function change_parameter(orig_params_matrix, parameter_name, value, position)
    result = deepcopy(orig_params_matrix)
    result[position...] = merge(orig_params_matrix[position...], [parameter_name => value])
    return result
end

function playout_parameter_injection(init, params_matrix, parameter, value, change_position, obs_positions; show_gif = false, stop_t = 10, dir_path = "")
    new_p_matrix = change_parameter(params_matrix, parameter, value, change_position)
    new_result = spatial_simulation(init, new_p_matrix; stop_t = stop_t)
    visualize_simulation_results(new_result, obs_positions; show_gif = show_gif, dir_path = dir_path, comment = " with gamma = "*string(value)*" at position ["*string(change_position.x)*", "*string(change_position.y)*"]")
end

function randomize_parameter(params_matrix, parameter, min_value, max_value)
    for x in 1:size(params_matrix, 1)
        for y in 1:size(params_matrix, 2)
            params_matrix[x,y] = merge(params_matrix[x,y], [parameter => rand(min_value: max_value)])
        end
    end
end

end