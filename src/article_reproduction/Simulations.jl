module Simulations

export exploratory_analysis
export simple_simulation
export draw_phase_diagram


using NamedTupleTools
import PyPlot


"""
Creates a phase diagram of a given functions func given the parameters params.
"""
function draw_phase_diagram(orig_params, funct; minval = 0, maxval = 2, steps = 200, plot_num = nothing, title = "")
    plot_num = plot_num == nothing ? rand(1:1000000) : plot_num
    fn1(s,i) = funct(s, i, translate_nondimensional_scheme(orig_params), orig_params)  # we have the original parameters here to translate the nondim result back to the original result

    S = repeat(range(minval,stop=maxval,length=steps)',steps)
    I = repeat(range(minval,stop=maxval,length=steps),1,steps)

    DXY = fn1.(S, I)

    fig = PyPlot.figure("pyplot_streamplot_$plot_num",figsize=(5,5))
    PyPlot.subplot()
    PyPlot.streamplot(S,I,map(first,DXY),map(x->x[2], DXY))
    
    PyPlot.suptitle(title)
    
    PyPlot.xlabel("Susceptible")
    PyPlot.ylabel("Infested")
    return(fig)
end

"""
Saves figuer fig to the given path.
"""
function print_dia(fig, path)
    fig.savefig(path*".png")
end

"""
Creates non-dimensional parameters for the non-dimensionalized
    equations using the non-dimensionalization scheme in Table 3.
"""
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
    mig_out = orig_params.mig_out, 
    t = orig_params.t
    
    )   
    return return_params

end


"""
    exporatory_analysis(func, params, name, grid)

Runs simulation of given function func given the parameters params modified by cartesian product
    of parameters specified in grid and saves the phase diagram to file in the name path
    with a name describing the modified parameters.
"""
function exploratory_analysis(func, orig_params, filename, grid; results = [], save_to_files = false, title = "")

    pair = grid[1]
    key_name = pair[1]
    values = pair[2]

    for value in values
        new_filename = build_filename(filename, key_name, value)
        new_title = build_title(title, key_name, value)

        # orig = before non-dimensionalization
        new_orig_params = modify_params(orig_params, key_name, value)
        if length(grid) > 1
            grid_copy = copy(grid)
            popfirst!(grid_copy)
            
            exploratory_analysis(func, new_orig_params, new_filename, grid_copy, results = results, title = new_title, save_to_files = save_to_files) # recursion
            
        else
            # translation nondim

            fig = draw_phase_diagram(new_orig_params, func, title = new_title)

            if save_to_files
                print_dia(fig, new_filename)           
            end

            push!(results, (fig = fig, filename =new_filename))
        end
    end

    return results

end

function build_filename(orig_name, key_name, value)
    return orig_name * string(key_name) *"="* replace(string(value), '.' => '_') * ";"
end

function build_title(orig_title, key_name, value)
    if orig_title == ""
        return orig_title * string(key_name) *" = " *string(value)
    end
    return orig_title *", "* string(key_name) *" = " *string(value)
end

"""
Changes and returns one parameter in the parameters params.
Changes parameter named key_value to value.
Creates and returns a new NamedTuple instance, does not change the original instance. 
"""
function modify_params(params, key_name, value)
    return(merge(params, [key_name=>value]))    
end

"""
Runs a simple simulation from an initial condition and returns the playout as two vectors
"""
function simple_simulation(orig_params, init, functions; delta_t = 0.1, stop_t = 50)
    Ss = []
    Is = []

    S = init.S
    I = init.I

    for _ in 0:delta_t:stop_t
        push!(Ss, S)
        push!(Is, I)

        S, I = 
            S + functions.dS(S, I, translate_nondimensional_scheme(orig_params)) * delta_t,
            I + functions.dI(S, I, translate_nondimensional_scheme(orig_params)) * delta_t
    end
    return Ss, Is
end

end