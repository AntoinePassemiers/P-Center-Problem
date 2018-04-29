# -*- coding: utf-8 -*-
# p3.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using JuMP


function formulation_3(p::Int64,
                       solver::Any,
                       rho::Array{Int64},
                       a::Array{Int64})
    model = Model(solver=solver)
    N, M, T = size(a)

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


function create_rph(p::Int64,
                    solver::Any,
                    rho::Array{Int64},
                    a::Array{Int64},
                    h::Int64)
    model = Model(solver=solver)

    N, M, T = size(a)

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


function solve_p3(d::Array{Array{Int64}}, p::Int64, solver::Any)
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

    model, y = formulation_3(p, solver, rho, a)
    status = solve(model)
    return model, y, status
end


function solve_p3_with_BINARY(d::Array{Array{Int64}}, p::Int64, solver::Any)
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
    while max_h - min_h > 1
        mid::Int64 = convert(Int64, floor((min_h + max_h) / 2))
        RPh = create_rph(p, solver, rho, a, mid)
        status = solve(RPh)
        println(status)
        if status == :Infeasible
            min_h = mid
        else
            max_h = mid
            LB = rho[mid]
        end
    end

    rho = [x for x in Iterators.filter(x->x>=LB, rho)]
    T = length(rho)
    a = Array{Int64}(N, M, T)
    for i=1:N, j=1:M, k=1:T
        a[i, j, k] = (d[i][j] <= rho[k])
    end

    model, y = formulation_3(p, solver, rho, a)
    status = solve(model)
    return model, y, status
end



function solve_p3_with_DB3(d::Array{Array{Int64}}, p::Int64, solver::Any)
    # Flatten the distance matrix, keep only radii comprised between
    # UB and LB, and sort the values
    rho::Array{Int64} = unique(Iterators.flatten(d))
    sort!(rho)

    N::Int64 = length(d)    # Number of nodes
    M::Int64 = length(d[1]) # Number of possible centers
    T::Int64 = length(rho)  # Number of radius values

    a::Array{Int64} = Array{Int64}(N, M, T)
    for i = 1:N, j = 1:M, k = 1:T
        a[i, j, k] = (d[i][j] <= rho[k])
    end

    _min::Int64, _max::Int64 = 1, T
    while _max - _min >= 1
        # DB3
        _a::Int64 = _min + convert(Int64, floor((_max - _min) / 3))
        _b::Int64 = _min + 2 * convert(Int64, floor((_max - _min) / 3))
        println(size(a[:, :, [_a, _b]]))
        model, _ = formulation_3(p, solver, rho[[_a, _b]], a[:, :, [_a, _b]])
        status = solve(model)
        println(status)
        if status == :Infeasible
            _min = _max
            break
        else
            obj = convert(Int64, round(getobjectivevalue(model)))
            if obj == rho[_a]
                _max = _a
            elseif obj == rho[_b]
                _min = _a + 1
                _max = _b
            else
                _min = _b + 1
            end
        end
    end

    model, y = formulation_3(p, solver, rho[_min:_max], a[:, :, _min:_max])
    status = solve(model)
    return model, y, status
end