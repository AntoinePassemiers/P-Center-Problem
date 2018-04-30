# P-Center Problem (PCP) solver

## Dependencies

```julia
julia> Pkg.add("ArgParse")
julia> Pkg.add("Cbc")
julia> Pkg.add("GLPKMathProgInterface")
julia> Pkg.add("JuMP")
```

## How to run solve.jl

General command syntax:

```sh
$ julia solve.jl <path-to-instance> <form> [--solver SOLVER]
```

Argument "form" must be either "p1", "p3", "p3-binary" or "p3-db3".
Argument "solver" can be either "cbc" or "glpk" (default value is "cbc").

The script will create a folder "results" in current directory
and write the results for given formulation, instance and solver
in a txt file.

### Run formulation P1 the fastest way

```sh
$ julia solve.jl path-to-instance p1
```

### Run formulation P3 the fastest way

```sh
$ julia solve.jl path-to-instance p3-db3
```