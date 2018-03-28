# -*- coding: utf-8 -*-
# p1.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using JuMP
using Cbc


function solve_p1(d::Array{Array{Int64}}, p::Int64)
    N = length(distances)
    m = Model(solver=CbcSolver())

    # Define variables
    @variables m begin
        x[1:N,1:N], Bin
        y[1:N] >= 0, Bin
        z >= 0, Int
    end

    # Set objective
    @objective(m, Min, z)

    # Define constraints
    for i in 1:N
        # Constraint 2: sum_j (d_ij * x_ij) <= z  for each i=1...N
        @constraint(m, sum(d[i][j]*x[i, j] for j=1:N) <= z)

        # Constraint 3: sum_j (x_ij) == 1  for each i=1...N
        @constraint(m, sum(x[i, j] for j=1:N) == 1)

        # Constraint 4: x_ij <= y_j  for each i,j=1...N
        for j in 1:N
            @constraint(m, x[i, j] <= y[j])
        end
    end
    # Constraint 5: sum_j (y_j) <= p
    @constraint(m, sum(y) <= p)

    # Solve p-center problem
    status = solve(m)
    println("Status: ", status)

    z_value = getvalue(z)
    println(z_value)
    println(getvalue(x))

    for i in 1:N
        for j in 1:N

        end
    end

    return z_value
end
