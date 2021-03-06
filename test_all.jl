# -*- coding: utf-8 -*-
# test_all.jl: Run solver on all instances and save execution times
# authors : Antoine Passemiers, Cedric Simar

include("solve.jl")


implementations = ["p1", "p3-binary", "p3-db3"]
solvers = ["cbc", "glpk"]
EASY_ONLY = false


"""
    get_instance_filepaths(folder)

Return the paths of all .dat files located in folder.
Directory exploration is non-recursive.

# Arguments
- `folder::AbstractString`: path to the directory to explore.
"""
function get_instance_filepaths(folder::AbstractString)
    filenames = filter(x->contains(x, ".dat"), readdir(folder))
    return [joinpath(folder, filename) for filename in filenames]
end


"""
    save_Nx2matrix(mat, filenames, dest)

Save table of results at given location.

# Arguments
- `mat::Array{Float64}`: table of results (times or objective values).
- `filepaths::Any`: sequence of file paths. Its length must be the number
   of rows in mat, since each row represents the results for an instance file.
- `dest::AbstractString`: destination path
"""
function save_Nx2matrix(mat::Array{Float64},
                        filepaths::Any,
                        dest::AbstractString)
    nrows, ncols = size(mat)
    open(dest, "w") do f
        write(f, string(implementations), "\r\n")
        for i = 1:nrows
            name = basename(filepaths[i])
            write(f, "$name    ")
            for j = 1:length(implementations)
                write(f, string(mat[i, j]), "    ")
            end
            write(f, "\r\n")
        end
    end
end


function main()
    # Get names of easy instance files and hard instance files
    easy_files = get_instance_filepaths("../instances/easy")
    if EASY_ONLY
        files = easy_files
    else
        hard_files = get_instance_filepaths("../instances/hard")
        files = [easy_files; hard_files]
    end

    # Allocate arrays for storing results
    n_instances = length(files)
    exectimes = Array{Float64}(n_instances, length(implementations))
    obj = Array{Float64}(n_instances, length(implementations))

    # Create results folder if not exists
    if isdir("../results") == false
        mkdir("../results")
    end

    # Store results for each solver and instance
    for solver in solvers
        for i = 1:n_instances, k = 1:length(implementations)
            parameters = Dict{String, Any}(
                "filepath"=>files[i],
                "form"=>implementations[k],
                "solver"=>solver)
            model, y, exectime = solve_p_center(parameters)
            exectimes[i, k] = exectime
            obj[i, k] = round(getobjectivevalue(model))
        end
        save_Nx2matrix(exectimes, files, "../results/times_$solver.txt")
        save_Nx2matrix(obj, files, "../results/obj_$solver.txt")
    end
end


main()