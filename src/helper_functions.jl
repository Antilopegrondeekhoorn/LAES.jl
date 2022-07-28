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

function isentropic_T(state_in,p_out)
    # a guess of the isentropic temperature at higher/lower pressure
    if ustrip(p_out) > state_in.p
        T_out_is_guess = (state_in.T+100)K
    else
        T_out_is_guess = (state_in.T-100)K
    end
    iterations = 0
    prev_steps = []
    while true
        iterations +=1
        s_is_guess = CoolProp.PropsSI("S","T",T_out_is_guess,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.frac_N2)]&Oxygen[$(1-state_in.frac_N2)]")  #composition doesn't change
        #change step when closer to solution
        diff = abs(ustrip(s_is_guess)-state_in.s) 
        if diff > 200
            step = 2K
        elseif 100 < diff < 200
            step = 1K
        elseif 25 < diff < 100
            step = 0.5K
        elseif 10 < diff < 25
            step = 0.1K
        elseif 2 < diff <10
            step = 0.05K
        else
            step = 0.01K
        end
        #change the temperature guess
        if diff < 0.05
            global T_out_is = T_out_is_guess
            break
        elseif ustrip(s_is_guess) < state_in.s
            T_out_is_guess += step
            push!(prev_steps,+1)
        elseif ustrip(s_is_guess) > state_in.s
            T_out_is_guess -= step
            push!(prev_steps,-1)
        end
        #prevent being in an infinite loop
        if iterations%4 == 0
            if sum(prev_steps) ==  0
                @warn "Infinite loop: adjust the steps"
                break
            end
            prev_steps = []
        end   
    end
    return T_out_is
end

function T_real(T_isentropic,state_in,p_out,h_out_real)
    T_out_real_guess = T_isentropic
    iterations = 0
    prev_steps = []
    while true
        iterations +=1
        h_out_real_guess = CoolProp.PropsSI("H","T",T_out_real_guess,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.frac_N2)]&Oxygen[$(1-state_in.frac_N2)]")  #composition doesn't change
        #change step when closer to solution
        diff = abs(ustrip(h_out_real_guess)-h_out_real) 
        if diff > 5000
            step = 2K
        elseif 1000 < diff < 5000
            step = 1K
        elseif 500 < diff < 1000
            step = 0.5K
        elseif 200 < diff < 500
            step = 0.1K
        elseif 100 < diff < 200
            step = 0.05K
        elseif 50 < diff < 100
            step = 0.01K
        else
            step = 0.001K
        end
        #change the temperature guess
        if diff < 1
            global T_out_real = T_out_real_guess
            break
        elseif ustrip(h_out_real_guess) < h_out_real
            T_out_real_guess += step
            push!(prev_steps,+1)
        elseif ustrip(h_out_real_guess) > h_out_real
            T_out_real_guess -= step
            push!(prev_steps,-1)
        end
        #prevent being in an infinite loop
        if iterations%4 == 0
            if sum(prev_steps) ==  0
                @warn "Infinite loop: adjust the steps"
                break
            end
            prev_steps = []
        end   
    end
    return T_out_real
end