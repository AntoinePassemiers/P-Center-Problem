# -*- coding: utf-8 -*-
# p3.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using JuMP


function create_rph(p::Int64, 
                    solver::Any,
                    rho::Array{Int64},
                    a::Array{Int64},
                    h::Int64)
    model = Model(solver=solver)

    N::Int64 = length(d)    # Number of nodes
    M::Int64 = length(d[1]) # Number of possible centers
    T::Int64 = length(rho)

    z::Array{Int64} = zeros(Int64, T)
    z[h] = 1

    # Define variables
    @variables model begin
        y[1:M] >= 0
    end

    # Define objective function
    @objective(model, Min, sum(rho[k]*z[k] for k=1:T))

    # Define constraints
    for i = 1:N
        @constraint(model, sum(a[i, j, h]*y[j] for j=1:M) >= 1)
    end
    # Constraint 16
    #    sum_j y[j] <= p
    #    There are at most p selected centers
    @constraint(model, sum(y[j] for j=1:M) <= p)

    return model
end


function create_p3(d::Array{Array{Int64}}, p::Int64, solver::Any)
    model = Model(solver=solver)

    # Flatten the distance matrix, keep only radii comprised between
    # UB and LB, and sort the values
    rho::Array{Int64} = unique(Iterators.flatten(d))
    sort!(rho)

    N::Int64 = length(d)    # Number of nodes
    M::Int64 = length(d[1]) # Number of possible centers
    T::Int64 = length(rho)

    a::Array{Int64} = Array{Int64}(N, M, T)
    for i = 1:N, j = 1:M, k = 1:T
        a[i, j, k] = (d[i][j] <= rho[k])
    end

    # Apply BINARY algorithm to find a good lower bound
    min_h::Int64 = 1
    max_h::Int64 = T
    LB::Int64 = typemax(Int64)
    while max_h - min_h >= 1
        mid::Int64 = convert(Int64, floor((min_h + max_h) / 2))
        RPh = create_rph(p, solver, rho, a, mid)
        status = solve(RPh)
        println(status)
        if status == "Infeasible"
            min_h = mid
        else
            max_h = mid
            LB = rho[mid]
        end
        println(mid)
        println(LB)
    end


    rho = [x for x in Iterators.filter(x->x>=LB, rho)]
    T = length(rho)
    a = Array{Int64}(N, M, T)
    for i = 1:N, j = 1:M, k = 1:T
        a[i, j, k] = (d[i][j] <= rho[k])
    end

    println(rho)

    # Define variables
    @variables model begin
        y[1:M], Bin
        z[1:T], Bin
    end

    # Define objective function
    @objective(model, Min, sum(rho[k]*z[k] for k=1:T))

    # Define constraints
    for i = 1:N, k = 1:T
        @constraint(model, sum(a[i, j, k]*y[j] for j=1:M) >= z[k])
    end
    @constraint(model, sum(y[j] for j=1:M) <= p)
    @constraint(model, sum(z[k] for k=1:T) == 1)

    return model, y
end