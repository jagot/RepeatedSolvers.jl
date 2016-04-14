# RepeatedSolvers

Simple layer for solvers that are used to solve the same linear system
repeatedly, but with varying right-hand sides. Useful in
time-stepping:

    Ax⁽ⁱ⁺¹⁾ = x⁽ⁱ⁾

After construction, they all follow the same interface

```julia
solver = DirectSolver(A) # or JacobiSolver, &c
solve!(solver, x)
```

[Example usage](Example.ipynb)

[![Build Status](https://travis-ci.org/jagot/RepeatedSolvers.jl.svg?branch=master)](https://travis-ci.org/jagot/RepeatedSolvers.jl)
