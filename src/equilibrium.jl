
using JSON
using Optim

include("equations.jl")

"""
The parameter, which needs to be predicted.
"""
@enum Predict TEMPERATURE PRESSURE

"""
    EquilibriumPoint(eq, predict, phase; max_P=nothing, max_T=nothing,
                     tol=nothing, incr=nothing, incr_min=nothing,
                     incr_den=nothing)

Contains the data that necessary for predicting an equilibrium point using an
equation of state.
"""
type EquilibriumPoint

    eq::Equation
    predict::Predict
    phase::Phase
    max_P::Float64
    max_T::Float64
    tol::Float64
    incr::Float64
    incr_min::Float64
    incr_den::Float64

    function EquilibriumPoint(eq, predict, phase; max_P=nothing, max_T=nothing,
                              tol=nothing, incr=nothing, incr_min=nothing,
                              incr_den=nothing)

        max_P == nothing && (max_P = max(eq.mix.Pcr...) * 1.1)
        max_T == nothing && (max_T = max(eq.mix.Tcr...) * 1.3)
        tol == nothing && (tol = 0.001)
        if predict == PRESSURE
            incr == nothing && (incr = 101325.0)
            incr_min == nothing && (incr_min = 500.0)
            incr_den == nothing && (incr_den = 200.0)
        else
            incr == nothing && (incr = 10.0)
            incr_min == nothing && (incr_min = 0.05)
            incr_den == nothing && (incr_den = 2.0)
        end

        return new(eq, predict, phase, max_P, max_T, tol, incr,
                   incr_min, incr_den)
    end
end

"""
    CriticalPoint(ep; T=nothing, P=nothing, incr=nothing, incr_min=nothing,
                  incr_den=nothing, corr_coef=nothing)

Contains the data that necessary for predicting an critical point using an
equation of state.
"""
type CriticalPoint

    ep::EquilibriumPoint
    T::Float64
    P::Float64
    incr::Float64
    incr_min::Float64
    incr_den::Float64
    corr_coef::Float64

    function CriticalPoint(ep; T=nothing, P=nothing, incr=nothing,
                           incr_min=nothing, incr_den=nothing,
                           corr_coef=nothing)
        T == nothing && (T = sum(ep.eq.mix.Tb .* ep.eq.mix.x))
        P == nothing && (P = 101325.0)
        if ep.predict == PRESSURE
            incr == nothing && (incr = 10.0)
            incr_min == nothing && (incr_min = 0.05)
            incr_den == nothing && (incr_den = 14.14)
            corr_coef == nothing && (corr_coef = 0.9)
        else
            incr == nothing && (incr = 202650.0)
            incr_min == nothing && (incr_min = 3000.0)
            incr_den == nothing && (incr_den = 8.21)
            corr_coef == nothing && (corr_coef = 0.99)
        end

        return new(ep, T, P, incr, incr_min, incr_den, corr_coef)

    end
end

################################################################################
## Equilibrium Point
################################################################################

"""
    solve(ep::EquilibriumPoint, T::Float64, P::Float64)

Solves an equation of state and returns an equilibrium points for the given
properties.
"""
function solve(eps_::Array{EquilibriumPoint, 1}, T::Float64, P::Float64)
	function app(a::Array{EquilibriumData, 1}, b::Array{EquilibriumData, 1})
        append!(a, b)
        return a
    end
    eds = @parallel app for ep ∈ eps_
        [solve(ep, T, P)]
    end
    return eds
end

"""
    solve(ep::EquilibriumPoint, T::Float64, P::Float64)

Solves an equation of state and returns an equilibrium point for the given
properties.
"""
function solve(ep::EquilibriumPoint, T::Float64, P::Float64)
    #incr = ep.incr
    ed = predict!(deepcopy(ep), T, P)
    #ep.incr = incr
    return ed
end

######################
## TEST VERSION
######################

function solve2(ep::EquilibriumPoint, T::Float64, P::Float64)
    x = copy(ep.eq.mix.x)
    y = copy(ep.eq.mix.y)
    ed = predict2!(ep, T, P)
    ep.eq.mix.x = x
    ep.eq.mix.y = y
    return ed
end

function find_computable(ep::EquilibriumPoint, T::Float64, P::Float64, crit::Float64, slv, val::Float64)
    
    incr = ep.incr
    incr_den = ep.incr_den
    incr_min = ep.incr_min
    Rval = val
    while incr > incr_min
        ed = slv(val)
        if ed == nothing
            if val > crit
                incr /= incr_den
                val = Rval
            end
            val += incr
        else
            return ed, val
        end
    end
    return nothing, 0
    
end

function predict2!(ep::EquilibriumPoint, T::Float64, P::Float64)

    if ep.predict == PRESSURE
        crit = ep.max_P
        slv = (x) -> solve(ep.eq, T, x)
        adj_s = (x) -> x
        val = P
    else
        crit = ep.max_T
        slv = (x) -> solve(ep.eq, x, P)
        adj_s = (x) -> 1.1 - x / 10.0
        val = T
    end
    
    ed, val = find_computable(ep, T, P, crit, slv, val)
    if ed == nothing return nothing end
    is_vapor = ep.phase == VAPOR
    rmf = is_vapor ? rmfV : rmfL
    p_val = val
    s = 0
    ns = 0
    while abs(1 - s) > ep.tol
        ed = slv(val)
        if ed == nothing
            x = ns - 1
            if abs(x) < 1e-5 return nothing end
            ns = 1 + x/2.5
            val = is_vapor ? p_val/ns : p_val*ns
        else
            p_val = val
            s = adj_s(rmf(ep, ed))
            ns = s
            val = is_vapor ? val/s : val*s
        end
    end
    return ed
    
end

function rmfV(ep, ed)
    k = ed.ϕv ./ ed.ϕl
    kmf = k .* ed.y
    s = sum(kmf)
    nx = kmf/s
    ed.x = nx
    return s
end

function rmfL(ep, ed)
    k = ed.ϕl ./ ed.ϕv
    kmf = k .* ed.x
    s = sum(kmf)
    ny = kmf/s
    ed.y = ny
    return s 
end

######################
## TEST VERSION
######################

"""
    predict!(ep::EquilibriumPoint, T::Float64, P::Float64)

Solves an equation of state and returns an equilibrium point for the given
properties. But it chenges the increment.
"""
function predict!(ep::EquilibriumPoint, T::Float64, P::Float64)

    s = 0
    if ep.predict == PRESSURE
        crit = ep.max_P
        slv = (x) -> solve(ep.eq, T, x)
        adj_s = (x) -> x
        val = P
    else
        crit = ep.max_T
        slv = (x) -> solve(ep.eq, x, P)
        adj_s = (x) -> 1.1 - x / 10.0
        val = T
    end
    x::Array{Float64, 1} = copy(ep.eq.mix.x)
    y::Array{Float64, 1} = copy(ep.eq.mix.y)
    Rval = val
    while true
        ed = slv(val)
        if ed == nothing
            if val > crit
                ep.incr /= ep.incr_den
                val = Rval
            end
            val += ep.incr
        else
            if ep.phase == VAPOR
                s = adj_s(update_x!(ep, ed))
                val /= s
            else
                s = adj_s(update_y!(ep, ed))
                val *= s
            end
        end
        ep.incr < ep.incr_min && return nothing
		if abs(1 - s) < ep.tol return ed end
    end
end

"""
    update_x!(ep, ed)

Updates liquid fraction in the given mixture and returns the S coefficient for
the further computations.
"""
function update_x!(ep, ed)
    k = ed.ϕv ./ ed.ϕl
    ky = k .* ed.y
    s = sum(ky)
    ep.eq.mix.x = ky / s
    return s
end

"""
    update_y!(ep, ed)

Updates vapor fraction in the given mixture and returns the S coefficient for
the further computations.
"""
function update_y!(ep, ed)
    k = ed.ϕl ./ ed.ϕv
    kx = k .* ed.x
    s = sum(kx)
    ep.eq.mix.y = kx / s
    return s
end

################################################################################
## Critical Point
################################################################################

#########################
## TEST VERSION
#########################

function solve2(cps::Array{CriticalPoint, 1})
    function app(a::Array{Array{EquilibriumData, 1}, 1},
                b::Array{Array{EquilibriumData, 1}, 1})
        append!(a, b)
        return a
    end
    edss = @parallel app for cp ∈ cps
        [solve2(cp)]
    end
    return edss
end

function solve2(cp::CriticalPoint)
    return predict2!(deepcopy(cp))
end

function predict2!(cp::CriticalPoint)
    
    eds::Array{EquilibriumData, 1}, ed = [], 1
    if cp.ep.predict == PRESSURE
        upd1! = (cp::CriticalPoint) -> cp.T -= cp.incr
        upd2! = (cp::CriticalPoint) -> cp.T += cp.incr
        upd3! = (cp::CriticalPoint) -> cp.P *= cp.corr_coef
    else
        upd1! = (cp::CriticalPoint) -> cp.P -= cp.incr
        upd2! = (cp::CriticalPoint) -> cp.P += cp.incr
        upd3! = (cp::CriticalPoint) -> cp.T *= cp.corr_coef
    end

    while ed != nothing
        ed = solve2(cp.ep, cp.T, cp.P)
        if ed != nothing
            println("T: $(ed.T), P: $(ed.P), X: $(cp.ep.eq.mix.x[1])")
            cp.T, cp.P = ed.T, ed.P
            push!(eds, ed)
            upd3!(cp)
        elseif cp.incr / cp.incr_den > cp.incr_min
            ed = 1
            upd1!(cp)
            cp.incr /= cp.incr_den
        end
        upd2!(cp)
    end

    return eds
end

#########################
## TEST VERSION
#########################

"""
    solve(cp::CriticalPoint)

Solves an equation of state and returns an array of equilibrium points, where
last of them a critical, for the given properties.
"""
function solve(cp::CriticalPoint)
    return predict!(deepcopy(cp))
end

"""
    predict!(cp::CriticalPoint)

Solves an equation of state and returns an an array of equilibrium points, where
last of them a critical, for the given properties. But is changes the original
properties.
"""
function predict!(cp::CriticalPoint)

    eds::Array{EquilibriumData, 1}, ed = [], 1
    if cp.ep.predict == PRESSURE
        upd1! = (cp::CriticalPoint) -> cp.T -= cp.incr
        upd2! = (cp::CriticalPoint) -> cp.T += cp.incr
        upd3! = (cp::CriticalPoint) -> cp.P *= cp.corr_coef
    else
        upd1! = (cp::CriticalPoint) -> cp.P -= cp.incr
        upd2! = (cp::CriticalPoint) -> cp.P += cp.incr
        upd3! = (cp::CriticalPoint) -> cp.T *= cp.corr_coef
    end

    while ed != nothing
        ed = solve(cp.ep, cp.T, cp.P)
        if ed != nothing
            cp.T, cp.P = ed.T, ed.P
            push!(eds, ed)
            upd3!(cp)
        elseif cp.incr / cp.incr_den > cp.incr_min
            ed = 1
            upd1!(cp)
            cp.incr /= cp.incr_den
        end
        upd2!(cp)
    end

    return eds
end

"""
    solve(cps::Array{CriticalPoint, 1})

Solves multiple critical points in parallel (if multiple threads is defined).
"""
function solve(cps::Array{CriticalPoint, 1})
    function app(a::Array{Array{EquilibriumData, 1}, 1},
                b::Array{Array{EquilibriumData, 1}, 1})
        append!(a, b)
        return a
    end
    edss = @parallel app for cp ∈ cps
        [solve(cp)]
    end
    return edss
end

"""
    fit(cps::Array{CriticalPoint, 1},
        k_bounds::Tuple{Float64,Float64}=(-0.2,0.2); Tcr_exp=nothing,
        Pcr_exp=nothing, Vcr_exp=nothing, w=[1.0,1.0,1.0], opt_set...)

Optimizes binary interaction parameter for given experimental critical data.
You may provide additional arguments to the optimization function (from Optim
using package) using keyword arguments.
"""
function fit(cps::Array{CriticalPoint, 1},
             k_bounds::Tuple{Float64,Float64}=(-0.2,0.2); Tcr_exp=nothing,
             Pcr_exp=nothing, Vcr_exp=nothing, w=[1.0,1.0,1.0], opt_set...)

    prms = [p for p ∈ [Tcr_exp, Pcr_exp, Vcr_exp] if p != nothing]
    if length(prms) == 0
        error("You have to specify at least one of the critical parametes.")
    end
    function to_opt(k::Float64)
        for cp ∈ cps
            cp.ep.eq.mix.k = [0.0 k; k 0.0]
        end
        edss = solve(cps)
        Tcr_calc::Array{Float64, 1} = [x[end].T for x ∈ edss]
        Pcr_calc::Array{Float64, 1} = [x[end].P for x ∈ edss]
        Vcr_calc::Array{Float64, 1} = [x[end].Vl for x ∈ edss]
        s::Float64 = 0.0
        for (i, (e, c)) ∈ enumerate(zip([Tcr_exp, Pcr_exp, Vcr_exp],
                                      [Tcr_calc, Pcr_calc, Vcr_calc]))
            if e != nothing
                s += sum(w[i] * ((e - c) ./ e).^2)
            end
        end
        return s
    end
    return optimize(to_opt, k_bounds[1], k_bounds[2]; show_trace=true,
                    opt_set...)
end

#########################
## TEST VERSION
#########################

function fit2(cps::Array{CriticalPoint, 1},
             k_bounds::Tuple{Float64,Float64}=(-0.2,0.2); Tcr_exp=nothing,
             Pcr_exp=nothing, Vcr_exp=nothing, w=[1.0,1.0,1.0], opt_set...)

    prms = [p for p ∈ [Tcr_exp, Pcr_exp, Vcr_exp] if p != nothing]
    if length(prms) == 0
        error("You have to specify at least one of the critical parametes.")
    end
    function to_opt(k::Float64)
        for cp ∈ cps
            cp.ep.eq.mix.k = [0.0 k; k 0.0]
        end
        edss = solve2(cps)
        Tcr_calc::Array{Float64, 1} = [x[end].T for x ∈ edss]
        Pcr_calc::Array{Float64, 1} = [x[end].P for x ∈ edss]
        Vcr_calc::Array{Float64, 1} = [x[end].Vl for x ∈ edss]
        s::Float64 = 0.0
        for (i, (e, c)) ∈ enumerate(zip([Tcr_exp, Pcr_exp, Vcr_exp],
                                      [Tcr_calc, Pcr_calc, Vcr_calc]))
            if e != nothing
                s += sum(w[i] * ((e - c) ./ e).^2)
            end
        end
        return s
    end
    return optimize(to_opt, k_bounds[1], k_bounds[2]; show_trace=true,
                    opt_set...)
end

#########################
## TEST VERSION
#########################

################################################################################
## PostProcessing Functions
################################################################################

"""
    dump_json(eds::Array{EquilibriumData, 1}, filename::AbstractString)

Saves the given array of equilibrium points in json format.
"""
function dump_json(eds::Array{EquilibriumData, 1}, filename::AbstractString)

    typedict(x) = Dict(fn=>getfield(x, fn) for fn ∈ fieldnames(x))
    to_json = [typedict(x) for x ∈ eds]
    open(filename, "w") do file
        txt = JSON.json(to_json)
        write(file, txt)
    end
end

"""
    dump_json(edss::Array{Array{EquilibriumData, 1}, 1},
                   filename::AbstractString)

Saves the given array of arrays of equilibrium points in json format.
The main use is for saving the product of multiple critical point valuations.
"""
function dump_json(edss::Array{Array{EquilibriumData, 1}, 1},
                   filename::AbstractString)

    typedict(x) = Dict(fn=>getfield(x, fn) for fn ∈ fieldnames(x))
    to_json = [[typedict(x) for x ∈ y] for y ∈ edss]
    open(filename, "w") do file
        txt = JSON.json(to_json)
        write(file, txt)
    end
end
