module Thermo

__precompile__()

export EquilibriumData, Component, Mixture
export SRK, PR, RK, VdW, Equation
export pressure, temperature, liquid, vapor
export CriticalPoint, EquilibriumPoint
export solve, dump_json, solve_multiple, fit

# include("equations.jl")
include("equilibrium.jl")

end
