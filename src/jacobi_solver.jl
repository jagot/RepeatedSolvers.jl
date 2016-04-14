# The Jacobi method for solving an equation system
# Ax = b
# is an iterative method, where the k:th approximation
# to the solution is calculated as
# x⁽ᵏ⁺¹⁾ = D⁻¹(b - Rx⁽ᵏ⁾),
# where D is the diagonal elements of A, and
# R = A - D.

type JacobiSolver <: Solver
    A
    Dinv
    R
    tol
    ω
    maxiter
    function JacobiSolver(A, tol, ω=2./3, maxiter = 10)
        D = Diagonal(diag(A))
        new(A, inv(D), A-D, tol, ω, maxiter)
    end
end

# This implementation assumes that x is reasonably "close" to b; this
# is the case if A = I + Ξ, with |Ξ| << 1.
function solve!(b::StridedVecOrMat, J::JacobiSolver)
    xp = ones(b)
    x = xp
    for k = 1:J.maxiter
        x = J.ω*(J.Dinv*(b - J.R*xp)) + (1-J.ω)*xp
        if norm(J.A*x-b) < J.tol
            break
        end
        xp = x
    end
    b[:] = x
end

export JacobiSolver, solve!
