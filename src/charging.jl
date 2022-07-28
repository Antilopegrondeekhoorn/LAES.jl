# Packages
using CoolProp
using Unitful: m, J, kPa, Pa, kg, K, W, °C, l, bar, MW, g, uconvert, °, ustrip, s, kW, hr,mol,kJ
using Plots
using DataFrames
using CurveFit

# State struct
mutable struct State
    p::Float64
    T::Float64
    h::Float64
    s::Float64
    mdot::Float64
    frac_N2::Float64
    function State(p,T,frac_N2,mdot)
        s = new()
        s.p = p
        s.T = T
        s.frac_N2 = frac_N2
        frac_O2 = 1-frac_N2
        s.h = CoolProp.PropsSI("H","T",T,"P|gas",p,"PR::Nitrogen[$frac_N2]&Oxygen[$frac_O2]")
        s.s = CoolProp.PropsSI("S","T",T,"P|gas",p,"PR::Nitrogen[$frac_N2]&Oxygen[$frac_O2]")
        s.mdot = mdot
        return s
    end
end

