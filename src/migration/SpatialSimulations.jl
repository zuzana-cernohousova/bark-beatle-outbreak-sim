module SpatialSimulations

export spatial_simulation
export translate_nondimensional_scheme
using ..MigrationSimFunctions

function model_step(matrix_S, matrix_I, params_matrix, k)
    mig_directions = [[0,1], [1, 0], [0, -1], [-1, 0]]

    emigrants = matrix_I
    immigrants = get_all_immigrant_numbers(emigrants, params_matrix, mig_directions)

    dS = dS_mig.(matrix_S, matrix_I, immigrants - emigrants, params_matrix)
    dI = dI_mig.(matrix_S, matrix_I, immigrants - emigrants, params_matrix)

    new_Ss, new_Is = deepcopy(matrix_S) + dS *k, deepcopy(matrix_I) + dI*k

    return new_Ss, new_Is
end

function translate_nondimensional_scheme(orig_params)
    # the non-dim parameters are on the left side, same parameter names are used to simplify the equations

    return_params = (
    # vital dynamics of trees
    g = orig_params.g/orig_params.d,
    K = orig_params.lambda*orig_params.K/orig_params.m,

    # tree infestation rate
    beta = orig_params.beta/orig_params.d,
    gamma = orig_params.r*orig_params.gamma/orig_params.e,
    beta_function = orig_params.beta_function,

    # aggregation
    n = orig_params.n,

    # migration
    mu = orig_params.lambda*orig_params.mu/(orig_params.e*orig_params.m),
    a = orig_params.a
    )   
    return return_params

end

function nondim_params_matrix(params_matrix)
    return translate_nondimensional_scheme.(params_matrix)
end

function spatial_simulation(init, orig_params_matrix; stop_t = 50, k = 0.1)
    progress_S = [init.S]
    progress_I = [init.I]
    
    index = 1
    for _ in 0:k:stop_t
        S, I = model_step(progress_S[index], progress_I[index], nondim_params_matrix(orig_params_matrix), k)
        
        push!(progress_S, S)
        push!(progress_I, I)
        index += 1
    end
    return progress_S, progress_I
end

function get_all_immigrant_numbers(emigrants, params_matrix, mig_directions)
    immigrants = deepcopy(emigrants)

    for x in 1:size(immigrants,1)
        for y in 1:size(immigrants,2)
            im = 0
            for dir in mig_directions
                im = im + emigrants[mod(x+dir[1]-1, size(immigrants, 1))+1 , mod(y+dir[2]-1, size(immigrants,2))+1]
            end

            immigrants[x, y] = im/length(mig_directions)
        end
    end
    return immigrants
end

end