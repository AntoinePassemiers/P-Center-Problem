# -*- coding: utf-8 -*-
# test_all.jl
# authors : Antoine Passemiers, Cedric Simar

include("solve.jl")


function get_instance_filepaths(folder)
    filenames = filter(x->contains(x, ".dat"), readdir(folder))
    return [joinpath(folder, filename) for filename in filenames]
end

formulations = ["p1", "p3"]
easy_files = get_instance_filepaths("instances/easy")
hard_files = get_instance_filepaths("instances/hard")
files = [easy_files; hard_files]
n_instances = length(files)

exectimes = Array{Float64}(n_instances, 2)
obj = Array{Int64}(n_instances, 2)

for i = 1:n_instances
    for k = 1:2
        model, y, exectime = solve_p_center(files[i], formulations[k])
        exectimes[i, k] = exectime
        obj[i, k] = round(getobjectivevalue(model))
    end
end

open("results/times.txt", "w") do f
    for i=1:n_instances
        t1 = exectimes[i, 1]
        t2 = exectimes[i, 2]
        write(f, "$t1    $t2 \r\n")
    end
end

open("results/obj.txt", "w") do f
    for i=1:n_instances
        o1 = obj[i, 1]
        o2 = obj[i, 2]
        write(f, "$o1    $o2 \r\n")
    end
end