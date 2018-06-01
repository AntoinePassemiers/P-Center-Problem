# -*- coding: utf-8 -*-
# p3.jl: Formulation P3 and algorithms for finding bounds for P3
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using JuMP

"""
    formulation_3(p, solver, rho, a)

Create a JuMP model and formulate P3. Return the model and an array y
where y[j] is a JuMP variable indicating whether center j is selected.

# Arguments
- `p::Int64`: the maximum number of selected centers.
- `solver::Any`: a reference to a JuMP-compatible solver (either Cbc or GLPK).
- `rho::Array{Int64}`: sorted unique values from distance matrix.
- `a::Array{Int64}`: array of dimensionality 3 where a[i, j, k] indicates
   whether distance d[i, j] is less or equal to radius value rho[k].
"""
function formulation_3(p::Int64, solver::Any, rho::Array{Int64}, a::Array{Int64})
    model = Model(solver=solver)
    N, M, K = size(a)

    # Define variables
    @variables model begin
        y[1:M], Bin
        z[1:K], Bin
    end

    # Define objective function
    @objective(model, Min, sum(rho[k]*z[k] for k=1:K))

    # Define constraints
    for i = 1:N, k = 1:K
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
    @constraint(model, sum(z[k] for k=1:K) == 1)
    return model, y
end


"""
    create_rph(p, solver, rho, a, h)

Create a JuMP model, formulate RP_h and return the model.

# Arguments
- `p::Int64`: the maximum number of selected centers.
- `solver::Any`: a reference to a JuMP-compatible solver (either Cbc or GLPK).
- `rho::Array{Int64}`: sorted unique values from distance matrix.
- `a::Array{Int64}`: array of dimensionality 3 where a[i, j, k] indicates
   whether distance d[i, j] is less or equal to radius value rho[k].
- `h::Int64`: Index of the variable to be set to 1.
"""
function create_rph(p::Int64,
                    solver::Any,
                    rho::Array{Int64},
                    a::Array{Int64},
                    h::Int64)
    model = Model(solver=solver)
    N, M, K = size(a)

    # Set all variables z[k] to 0 except for k == h
    z::Array{Int64} = zeros(Int64, K)
    z[h] = 1

    # Define variables
    # Linear relaxation on variable y
    @variables model begin
        y[1:M] >= 0
    end

    # Define objective function
    @objective(model, Min, sum(rho[k]*z[k] for k=1:K))

    # Define constraints
    for i = 1:N
        # Constraints group 23
        # Same constraints as in constraints group 15 except
        # that only k == h is considered
        @constraint(model, sum(a[i, j, h]*y[j] for j=1:M) >= 1)
    end
    # Constraint 24
    #    sum_j y[j] <= p
    #    There are at most p selected centers
    #    Same constraint as constraint 16
    @constraint(model, sum(y[j] for j=1:M) <= p)
    return model
end


"""
    get_p3_constants(d)

Compute all constants for formulation P3. More precisely,
compute matrix a, sort unique radius values (rho), 
and get the dimensions of the problem (N, M and T).
Everything is then returned.

# Arguments
- `d::Int64`: distance matrix of dimensionality 2.
"""
function get_p3_constants(d::Array{Array{Int64}})
    # Flatten the distance matrix, keep only radii comprised between
    # UB and LB, and sort the values
    rho::Array{Int64} = unique(Iterators.flatten(d))
    sort!(rho)

    N::Int64 = length(d)    # Number of nodes
    M::Int64 = length(d[1]) # Number of possible centers
    K::Int64 = length(rho)  # Number of radius values

    # Compute matrix a
    a::Array{Int64} = Array{Int64}(N, M, K)
    for i = 1:N, j = 1:M, k = 1:K
        a[i, j, k] = (d[i][j] <= rho[k])
    end
    return N, M, K, rho, a
end


"""
    solve_p3(d, p, solver)

Enumerate sorted unique radius values (rho), compute matrix a,
create a P3 formulation and solve it as an integer program.

Return the model, the array y of variables and the solve status.

# Arguments
- `d::Int64`: distance matrix of dimensionality 2.
- `p::Int64`: the maximum number of selected centers.
- `solver::Any`: a reference to a JuMP-compatible solver (either Cbc or GLPK).
"""
function solve_p3(d::Array{Array{Int64}}, p::Int64, solver::Any)
    # Compute matrix a and get unique radius values in ascending order
    N, M, K, rho, a = get_p3_constants(d)

    # Create model P3 and solve it
    model, y = formulation_3(p, solver, rho, a)
    status = solve(model)
    return model, y, status
end


"""
    solve_p3_with_BINARY(d, p, solver)

Enumerate sorted unique radius values (rho), compute matrix a,
create a P3 formulation and solve it as an integer program
using LB as lower bound. LB is found using BINARY algorithm.

Return the model, the array y of variables and the solve status.

# Arguments
- `d::Int64`: distance matrix of dimensionality 2.
- `p::Int64`: the maximum number of selected centers.
- `solver::Any`: a reference to a JuMP-compatible solver (either Cbc or GLPK).
"""
function solve_p3_with_BINARY(d::Array{Array{Int64}}, p::Int64, solver::Any)
    # Compute matrix a and get unique radius values in ascending order
    N, M, K, rho, a = get_p3_constants(d)

    pretty_print = (lb, minh, maxh) -> println(
        "    ", if (LB == typemax(Int64)) " " else "*" end,
        " min : ", rpad(rho[min_h], 5, " "), " <= ",
        "LB : ", rpad(if (LB == typemax(Int64)) "-" else LB end, 5, " "),
        " <= max : ", rpad(rho[max_h], 5, " "))

    # Apply BINARY algorithm to find a good lower bound
    min_h::Int64 = 1
    max_h::Int64 = K
    LB::Int64 = typemax(Int64) # LB = Infinity
    println("Processing best lower bound...")
    pretty_print(LB, min_h, max_h)
    while max_h - min_h > 1
        # Find mid by bisection search

        mid::Int64 = convert(Int64, floor((min_h + max_h) / 2))
        RPh = create_rph(p, solver, rho, a, mid)
        # Solve linear relaxation
        # Because we expect it to be infeasible at some iteration of BINARY,
        # warnings are disabled for estethic reasons.
        status = solve(RPh, suppress_warnings=true)
        if status == :Infeasible
            min_h = mid
        else
            max_h = mid
            LB = rho[mid]
        end
        pretty_print(LB, min_h, max_h)
    end

    # Remove all radius values below LB from the search space
    # and solve P3
    model, y = formulation_3(p, solver, rho[max_h:K], a[:, :, max_h:K])
    println("Solving original problem (LB = $LB)...")
    status = solve(model)
    return model, y, status
end


"""
    solve_p3_with_DB3(d, p, solver)

Enumerate sorted unique radius values (rho), compute matrix a,
create a P3 formulation and solve it using double bound (DB3) algorithm.

Return the model, the array y of variables and the solve status.

# Arguments
- `d::Int64`: distance matrix of dimensionality 2.
- `p::Int64`: the maximum number of selected centers.
- `solver::Any`: a reference to a JuMP-compatible solver (either Cbc or GLPK).
"""
function solve_p3_with_DB3(d::Array{Array{Int64}}, p::Int64, solver::Any)
    # Compute matrix a and get unique radius values in ascending order
    N, M, K, rho, a = get_p3_constants(d)

    pretty_print = (lb, ah, bh, ub) -> println(
        "    ",
        " LB : ", rpad(lb, 5, " "), " <= ",
        " a : ", rpad(ah, 5, " "), " < ",
        " b : ", rpad(bh, 5, " "), " <= ",
        "UB : ", rpad(ub, 5, " "))

    # Apply DB3 algorithm
    _min::Int64, _max::Int64 = 1, K
    println("Tightening bounds...")
    while _max - _min >= 1
        _a::Int64 = _min + convert(Int64, floor((_max - _min) / 3))
        _b::Int64 = _min + 2 * convert(Int64, floor((_max - _min) / 3))
        pretty_print(_min, _a, _b, _max)
        P_S, _ = formulation_3(p, solver, rho[[_a, _b]], a[:, :, [_a, _b]])
        status = solve(P_S, suppress_warnings=true)
        if status == :Infeasible
            _min = _b + 1
        else
            obj = convert(Int64, round(getobjectivevalue(P_S)))
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

    # Remove all radius values below _min and values above _max from 
    # the search space and solve P3
    model, y = formulation_3(p, solver, rho[_min:_max], a[:, :, _min:_max])
    println("Solving original problem (LB = $_min, UB = $_max)...")
    status = solve(model)
    return model, y, status
end