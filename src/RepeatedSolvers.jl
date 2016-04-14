module RepeatedSolvers

abstract Solver
include("direct_solver.jl")
include("jacobi_solver.jl")

end
