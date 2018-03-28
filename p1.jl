# -*- coding: utf-8 -*-
# p1.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using JuMP
using Cbc


function solve_p1(d::Array{Array{Int64}}, p::Int64)
    N = length(d)
    m = Model(solver=CbcSolver())

    # Define variables
    @variables m begin
        x[1:N,1:N], Bin
        y[1:N] >= 0, Bin
        z, Int
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
    println("Number of variables  : ", MathProgBase.numvar(m))
    println("Number of constraints: ", MathProgBase.numconstr(m))
    println("\nSolving problem...")
    status = solve(m)
    println("Status    : ", status)
    # println("Solve time: ", getsolvetime(m))

    y_values = getvalue(y)
    z_value = getvalue(z)
    println("Objective : ", z_value)
    return y_values, z_value
end
