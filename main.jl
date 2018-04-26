# -*- coding: utf-8 -*-
# main.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using ArgParse

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

function main()
    arg_settings = ArgParseSettings(description="P-Center Solver: parameters")

    @add_arg_table arg_settings begin
        "filepath"
            help = "Path to the instance file"
            required = true
        "form"
            help = "Formulation of the P-Center Problem"
            required = true
            default = "p1"
    end
    parsed_args = parse_args(ARGS, arg_settings)

    distances, p = load_instance(FILE_PATH)
    solver = CbcSolver()

    if parsed_args["form"] == "p1"
        model, y = create_p1(distances, p, solver)
    elseif parsed_args["form"] == "p3"
        model, y = create_p3(distances, p, solver)
    end

    # Solve p-center problem
    println("Number of variables  : ", MathProgBase.numvar(model))
    println("Number of constraints: ", MathProgBase.numconstr(model))
    println("\nSolving problem...")
    start_time = time()
    status = solve(model)
    exectime = time() - start_time
    println("Status    : $status")
    println("Solve time: $(@sprintf("%.3f", exectime)) s")

    obj = getobjectivevalue(model)
    println("Objective : ", obj)

    # Write solution
    open("out.txt", "w") do f
        write(f, "Value of the objective function: $obj \n\n")
        for i=1:length(y)
            if getvalue(y[i]) > 0
                write(f, "Center selected at area $i \n")
            end
        end
    end
end

main()