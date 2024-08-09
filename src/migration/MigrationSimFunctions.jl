module MigrationSimFunctions

export dI_mig
export dS_mig
export aggreg_beta_mig
export uniform_beta_mig

"""
    dI_mig(S, I, M, params)
    
Change of density of infested trees I per hectare.
M ... value of immigrants
"""
function dI_mig(S, I, M,  params)
    return(
        params.beta_function(S, I, M, params)*S
            # new infestation from spread and 
        - I # death from infestation
    )
end

"""
    dS(S, I, R, params)

Change of density of susceptible trees S per hectare.
"""
function dS_mig(S, I, M, params)
    return(
        G(S, I, params)         # new susceptible trees
        - params.beta_function(S, I, M, params)*S     # trees turning infected
        )
end


"""
    G(S,I, params)

Tree recruitment i.e. new susceptible trees
"""
function G(S,I, params) 
    return (params.g*(      # growth rate
            params.K        # maximum of trees / capacity of forest stand
            - S - I         # infected and susceptible trees do not contribute
            )
        )
end


"""
    aggreg_beta(S, I, params)

    Aggregated Pseudo steady infestation rate
"""
function aggreg_beta_mig(S, I, M, params)
    return(
        params.beta*(
            ((I + params.mu * migration_activation(M, params))^params.n)/(
                (I+ (params.mu * migration_activation(M, params)))^params.n 
                + params.gamma^params.n*(1 + S)^params.n
            )
        )
    )
end


"""
    uniform_beta_mig(S, I, M, params)

    Uniform Pseudo steady infestation rate
"""
function uniform_beta_mig(S, I, M, params)
    if (I + params.mu*migration_activation(M, params))/(1+S) < params.gamma
        return 0
    end
    return params.beta
end

function migration_activation(M, params)
    return 0.1 + (1.8/(1 + exp(-params.a*M)))
end


end