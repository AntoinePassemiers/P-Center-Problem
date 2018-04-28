# -*- coding: utf-8 -*-
# p3.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using JuMP


function create_p3(d::Array{Array{Int64}}, p::Int64, solver::Any)
    model = Model(solver=solver)

    N::Int64 = length(d)    # Number of nodes
    M::Int64 = length(d[1]) # Number of possible centers

    rho::Array{Int64} = unique([x for x in Iterators.flatten(d)])
    sort!(rho)
    T::Int64 = length(rho)

    a::Array{Int64} = Array{Int64}(N, M, T)
    for i = 1:N, j = 1:M, k = 1:T
        a[i, j, k] = (d[i][j] <= rho[k])
    end

    # Define variables
    @variables model begin
        y[1:M], Bin
        z[1:T], Bin
    end

    # Define objective function
    @objective(model, Min, sum(rho[k]*z[k] for k=1:T))

    # Define constraints
    for i = 1:N, k = 1:T
        # Constraints group 15
        #    sum_j a[i, j, k]*y[j] >= z[k]
        #    For each selected vertex, there is at least 
        #    one center that covers it within radius k
        @constraint(model, sum(a[i, j, k]*y[j] for j=1:M) >= z[k])
    end
    # Constraint 16
    #    sum_j y[j] <= p
    #    There are at most p selected centers
    @constraint(model, sum(y[j] for j=1:M) <= p)
    # Constraint 17
    #    sum_k z[k] == 1
    #    Exactly one distance D[k] is selected
    @constraint(model, sum(z[k] for k=1:T) == 1)

    return model, y

end