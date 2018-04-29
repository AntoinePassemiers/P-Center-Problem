# -*- coding: utf-8 -*-
# test_all.jl: Run solver on all instances and save execution times
# authors : Antoine Passemiers, Cedric Simar

include("solve.jl")


implementations = ["p1", "p3", "p3-binary", "p3-db3"]
solvers = ["cbc", "glpk"]
EASY_ONLY = false


function get_instance_filepaths(folder::AbstractString)
    filenames = filter(x->contains(x, ".dat"), readdir(folder))
    return [joinpath(folder, filename) for filename in filenames]
end

function save_Nx2matrix(mat::Array{Float64}, filepath::AbstractString)
    nrows, ncols = size(mat)
    open(filepath, "w") do f
        for i=1:nrows
            o1 = mat[i, 1]
            o2 = mat[i, 2]
            write(f, "$o1    $o2 \r\n")
        end
    end
end

easy_files = get_instance_filepaths("../instances/easy")
if EASY_ONLY
    files = easy_files
else
    hard_files = get_instance_filepaths("../instances/hard")
    files = [easy_files; hard_files]
end

n_instances = length(files)
exectimes = Array{Float64}(n_instances, 2)
obj = Array{Float64}(n_instances, 2)

# Create results folder if not exists
if isdir("../results") == false
    mkdir("../results")
end

for solver in solvers
    for i = 1:n_instances, k = 1:2
        parameters = Dict{String, Any}(
            "filepath"=>files[i],
            "form"=>implementations[k],
            "solver"=>solver)
        model, y, exectime = solve_p_center(parameters)
        exectimes[i, k] = exectime
        obj[i, k] = round(getobjectivevalue(model))
    end

    save_Nx2matrix(exectimes, "../results/times_$solver.txt")
    save_Nx2matrix(obj, "../results/obj_$solver.txt")
end