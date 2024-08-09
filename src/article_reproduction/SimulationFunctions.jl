module SimulationFunctions

export dI
export dS
export aggreg_beta
export uniform_beta
export simple_beta

"""
    dI(S, I, params)
    
Change of density of infested trees I per hectare.
"""
function dI(S, I, params)
    return(
        params.beta_function(S, I, params)*S
        - I
    )
end


"""
    dS(S, I, R, params)

Change of density of susceptible trees S per hectare.
"""
function dS(S, I, params)
    return (
        G(S, I, params)         # new susceptible trees
        - params.beta_function(S, I, params)*S     # trees turning infected
        )
end


"""
    G(S,I, params)

Tree recruitment i.e. new susceptible trees
"""
function G(S,I, params) 
    return (params.g*(      # growth rate
            params.K        # maximum of trees / capacity of forest stand
            - S - I         
            )
        )
end

"""
    aggreg_beta(S, I, params)

    Aggregated Pseudo steady infestation rate
"""
function aggreg_beta(S, I, params)
    return(
        params.beta*(
            ((I + params.mu)^params.n)/(
                (I+ params.mu)^params.n 
                + params.gamma^params.n*(1 + S)^params.n
            )
        )
    )
end


"""
    uniform_beta(S, I, params)

    Uniform Pseudo steady infestation rate
"""
function uniform_beta(S, I, params)
    if (I + params.mu)/(1+S) < params.gamma
        return 0
    end
    return params.beta
end

function simple_beta(S, I, params)
    return params.t *I
    
end

end