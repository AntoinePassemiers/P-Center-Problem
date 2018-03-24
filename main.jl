# -*- coding: utf-8 -*-
# main.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()    

using JuMP
using Cbc


INSTANCES_PATH = "instances"

FILE_PATH = joinpath(INSTANCES_PATH, "easy/instance10_1_1.dat")

distances = Array{Int64}[]
f = open(FILE_PATH)
lines = readlines(f)
header = split(lines[1], "\t")
n_sites = parse(Int64, header[1])
max_n_centers = parse(Int64, header[2])
for i = 1:n_sites
    elements = split(lines[1+i], "\t")
    row = map(i->parse(Int64, elements[i]), range(1, n_sites))
    push!(distances, row)
end
close(f)

println(distances)
println(max_n_centers)