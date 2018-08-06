module Thermo

__precompile__()

export EquilibriumData, Component, Mixture
export SRK, PR, RK, VdW, Equation
export PRESSURE, TEMPERATURE, LIQUID, VAPOR
export CriticalPoint, EquilibriumPoint
export solve, dump_json, fit

# include("equations.jl")
include("equilibrium.jl")

end
