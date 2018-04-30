# -*- coding: utf-8 -*-
# main.jl
# authors : Antoine Passemiers, Cedric Simar

__precompile__()

using ArgParse
using Cbc
using GLPKMathProgInterface

include("p1.jl")
include("p3.jl")


IMPLEMENTATIONS = ["p1", "p3", "p3-binary", "p3-db3"]

"""
    load_instance(file_path)

Parse a .dat file and return the distance matrix and
the maximum number of selected centers (p).

# Arguments
- `file_path::AbstractString`: Path to the instance file.
"""
function load_instance(file_path::AbstractString)
    distances = Array{Int64}[]
    f = open(file_path)
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


"""
    solve_p_center(parameters)

Given a solver, a formulation and an instance file,
create the corresponding problem, solve it and display
all useful information.

# Arguments
- `parameters::Dict{String, Any}`: parameters of PCP solver.
   Required keys are:
       solver::AbstractString: either "glpk" or "cbc".
       form::AbstractString: either "p1", "p3", "p3-binary" or "p3-db3".
       filepath::AbstractString: path to the instance file.
"""
function solve_p_center(parameters::Dict{String, Any})
    println("Loading instance from ", parameters["filepath"])
    distances, p = load_instance(parameters["filepath"])

    # Use chosen solver
    if parameters["solver"] == "cbc"
        solver = CbcSolver()
    elseif parameters["solver"] == "glpk"
        solver = GLPKSolverMIP()
    else
        error("Solver must be either cbc or glpk")
    end

    # Call appropriate algorithm and measure execution time
    start_time = time()
    if parameters["form"] == "p1"
        println("Using formulation p1")
        model, y, status = solve_p1(distances, p, solver)
    elseif parameters["form"] == "p3"
        println("Using formulation p3")
        model, y, status = solve_p3(distances, p, solver)
    elseif parameters["form"] == "p3-binary"
        println("Using formulation p3 and BINARY algorithm")
        model, y, status = solve_p3_with_BINARY(distances, p, solver)
    elseif parameters["form"] == "p3-db3"
        println("Using formulation p3 and DB3 algorithm")
        model, y, status = solve_p3_with_DB3(distances, p, solver)
    else
        error("Formulation must be either p1 or p3")
    end
    exectime = time() - start_time

    # Display useful information
    obj = convert(Int64, round(getobjectivevalue(model)))
    println("Number of variables  : ", MathProgBase.numvar(model))
    println("Number of constraints: ", MathProgBase.numconstr(model))
    println("Status    : $status")
    println("Solve time: $(@sprintf("%.3f", exectime)) s")
    println("Objective : $obj \n")

    # Create results folder if not exists
    if isdir("../results") == false
        mkdir("../results")
    end

    # Write solution
    results_path::AbstractString = replace(string(
        splitext(basename(parameters["filepath"]))[1], "_",
        parameters["form"], "_",
        parameters["solver"], ".txt"), "-", "_")
    open(joinpath("../results", results_path), "w") do f
        write(f, "Formulation: ", parameters["form"], "\r\n")
        write(f, "Solver: ", parameters["solver"], "\r\n")
        write(f, "Execution time: $exectime \r\n")
        write(f, "Value of the objective function: $obj \r\n\n")
        for i=1:length(y)
            if getvalue(y[i]) > 0
                write(f, "- Center selected at area $i \r\n")
            end
        end
    end
    return model, y, exectime
end


function main()
    # Parse command-line arguments
    arg_settings = ArgParseSettings(description="P-Center Solver: parameters")
    @add_arg_table arg_settings begin
        "filepath"
            help = "Path to the instance file"
            required = true
        "form"
            help = "Formulation of the P-Center Problem (p1 or p3)"
            required = true
            range_tester = (x->x in IMPLEMENTATIONS)
        "--solver"
            help = "Solver (either cbc or glpk)"
            required = false
            default = "cbc"
            range_tester = (x->x in ["cbc", "glpk"])
    end
    parsed_args = parse_args(ARGS, arg_settings)
    solve_p_center(parsed_args)
end


if length(ARGS) > 1
    main()
end