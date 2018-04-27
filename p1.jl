# -*- coding: utf-8 -*-
# p1.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using JuMP


function create_p1(d::Array{Array{Int64}}, p::Int64, solver)
    model = Model(solver=solver)

    N::Int64 = length(d)

    # Define variables
    @variables model begin
        x[1:N,1:N], Bin
        y[1:N] >= 0, Bin
        z, Int
    end

    # Set objective
    @objective(model, Min, z)

    # Define constraints
    for i in 1:N
        # Constraints group 2: sum_j (d_ij * x_ij) <= z  for each i=1...N
        #    The objective value is greater or equal to the
        #    maximum vertex-to-center distance.
        @constraint(model, sum(d[i][j]*x[i, j] for j=1:N) <= z)

        # Constraints group 3: sum_j (x_ij) == 1  for each i=1...N
        #    Each vertex is assigned to exactly one center.
        @constraint(model, sum(x[i, j] for j=1:N) == 1)

        # Constraints group 4: x_ij <= y_j  for each i,j=1...N
        #    A vertex can be assigned to v_j only if there is
        #    a center at v_j.
        for j in 1:N
            @constraint(model, x[i, j] <= y[j])
        end
    end
    # Constraint 5: sum_j (y_j) <= p
    #    There are at most p centers.
    @constraint(model, sum(y) <= p)

    return model, y
end
