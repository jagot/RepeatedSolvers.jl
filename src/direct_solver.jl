type DirectSolver <: Solver
    A
    function DirectSolver(A::AbstractMatrix)
        new(factorize(A))
    end
end

solve!(x, S::DirectSolver) =  A_ldiv_B!(S.A, x)

export DirectSolver, solve!
