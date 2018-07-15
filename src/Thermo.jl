module Thermo

__precompile__()

export EquilibriumData, Component, Mixture
export SRK, PR, RK, VdW
export pressure, temperature, liquid, vapor
export CriticalPoint, EquilibriumPoint
export solve, dump_json

include("equations.jl")
include("equilibrium.jl")

end
