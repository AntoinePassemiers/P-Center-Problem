# -*- coding: utf-8 -*-
# main.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

include("p1.jl")
include("p3.jl")


FILE_PATH = joinpath("instances", "easy/instance10_1_2.dat")


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
solver = CbcSolver()
model, x, y, z = create_p1(distances, p, solver)

# Solve p-center problem
println("Number of variables  : ", MathProgBase.numvar(model))
println("Number of constraints: ", MathProgBase.numconstr(model))
println("\nSolving problem...")
status = solve(model)
println("Status    : ", status)
# println("Solve time: ", getsolvetime(m))

println("Objective : ", getvalue(z))

# Write solution
open("out.txt", "w") do f
    write(f, "Value of the objective function: $z\n\n")
    for i=1:length(y)
        if getvalue(y[i]) > 0
            write(f, "Center at area $i\n")
        end
    end
end