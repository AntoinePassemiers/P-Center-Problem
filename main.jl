# -*- coding: utf-8 -*-
# main.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

include("p1.jl")
include("p3.jl")


FILE_PATH = joinpath("instances", "easy/instance10_1_1.dat")


function load_instance(file_path::AbstractString)
    distances = Array{Int64}[]
    f = open(FILE_PATH)
    lines = readlines(f)
    header = split(replace(lines[1], "\t", " "))
    n_sites = parse(Int64, header[1])
    max_n_centers = parse(Int64, header[2])
    for i = 1:n_sites
        elements = split(replace(lines[1+i], "\t", " "))
        row = map(i->parse(Int64, elements[i]), range(1, n_sites))
        push!(distances, row)
    end
    close(f)
    return distances, max_n_centers
end


distances, p = load_instance(FILE_PATH)
solve_p1(distances, p)
