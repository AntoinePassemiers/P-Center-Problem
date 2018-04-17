# -*- coding: utf-8 -*-
# p3.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using JuMP
using Cbc


function create_p3(d::Array{Array{Int64}}, p::Int64, solver)
    N = length(d)
    model = Model(solver=solver)


end