# Packages
using CoolProp
using Unitful: m, J, kPa, Pa, kg, K, W, °C, l, bar, MW, g, uconvert, °, ustrip, s, kW, hr,mol,kJ
using Plots
using DataFrames
using CurveFit

########################################################
# State struct
########################################################

mutable struct AirState
    fluid::String
    p::Float64
    T::Float64
    h::Float64
    s::Float64
    mdot::Float64
    frac_N2::Float64
    function AirState(p,T,mdot,frac_N2)
        s = new()
        s.fluid = "Air"
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
mutable struct OilState
    fluid::String
    p::Float64
    T::Float64
    h::Float64
    mdot::Float64
    function OilState(p,T,mdot)
        s = new()
        s.fluid = "Essotherm650"
        s.p = p
        s.T = T
        
        SpecificHeat = [0.0 1770; 20.0 1850; 40.0 1920; 60.0 1990; 80.0 2060; 100.0 2130; 150.0 2310; 200.0 2490; 250.0 2670; 300.0 2850; 320.0 2920]
        p = poly_fit(SpecificHeat[:,1] + ones(size(SpecificHeat)[1],1)*273.15,SpecificHeat[:,2], 2)
        cp(T) = p[1]+p[2]*T+p[3]*T^2
        h(T) = ((cp(T)+cp(273.15))/2)*(T-273.15)
        s.h = h(T)
  
        s.mdot = mdot

        return s
    end   
end

function State(fluid,p,T,mdot;frac_N2=nothing)
    if fluid == "Air"
        return AirState(p,T,mdot,frac_N2)
    elseif fluid == "Essotherm650"
        return OilState(p,T,mdot)
    else 
        @warn("Fluid not known")
    end
end
 
########################################################
# Other helper functions
########################################################

# Compressors and expanders

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

# Intercoolers

function T_Essotherm650(h)
    T = 288:0.001:800
    hs = []
    for i = 1:length(T)
        x = h_Essotherm650(T[i])
        push!(hs,x)
    end
    index = findfirst(x->x>h,hs)
    return T[index]
end

# Cryo-expander 

function interpolate_K(T,air_component,pressure_in_words)
    T = ustrip(T)
    if pressure_in_words == "1atm"
        Ks_N2 = [78 0.0798;80 0.3040;82 0.5282;84 0.7582;86 0.9766;88 1.2008;90 1.4249]
        Ks_O2 = [78 -1.13368;80 -1.1164;82 -0.8959;84 -0.6755;86 -0.4550;88 -0.2346;90 -0.0141]
    elseif pressure_in_words == "2atm"
        Ks_N2 = [84 0.7048;86 0.9030;88 1.1012;90 1.2994;92 1.4976;94 1.6958;96 1.8939]
        Ks_O2 = [84 -0.4573;86 -0.3016;88 -0.1459;90 0.0098;92 0.1655;94 0.3211;96 0.4768]
    elseif pressure_in_words == "5atm"
        Ks_N2 = [94 1.5503;96 1.7017;98 1.8531;100 2.0045;102 2.1559;104 2.3073;106 2.4588;108 2.6102]
        Ks_O2 = [94 0.6605;96 0.7877;98 0.9148;100 1.042;102 1.1692;104 1.2963;106 1.4235;108 1.5506]
    end
    
    if air_component == "Nitrogen"
        index = findfirst(x-> abs(x-T)<=2,Ks_N2[:,1])
        T1 = Ks_N2[index,1]
        T2 = Ks_N2[index+1,1]
        value1 = Ks_N2[index,2]
        value2 = Ks_N2[index+1,2]
    elseif air_component == "Oxygen"
        index = findfirst(x-> abs(x-T)<=2,Ks_O2[:,1])
        T1 = Ks_O2[index,1]
        T2 = Ks_O2[index+1,1]
        value1 = Ks_O2[index,2]
        value2 = Ks_O2[index+1,2]
    else 
        @warn("Component of the air not known. Air is considered a mixture of Nitrogen and Oxygen")
    end
    new_value = value1 + (T-T1)/(T2-T1)*(value2-value1)
    return new_value
end

function flash_calculation(p,frac_N2,yield,pressure_in_words) #p given with units
    Tbubble = CoolProp.PropsSI("T","P",p,"Q",0,"PR::Nitrogen[$(frac_N2)]&Oxygen[$(1-frac_N2)]")
    Tdew = CoolProp.PropsSI("T","P",p,"Q",1,"PR::Nitrogen[$(frac_N2)]&Oxygen[$(1-frac_N2)]")

    #Assume 80K (between Tbubble and Tdew)
    Tguess = ((round(ustrip(Tbubble))+round(ustrip(Tdew)))/2)*1.0K
    p = ustrip(p) #[Pa]
    while true
        global L_on_F = yield
        
        global K_nitrogen = exp(interpolate_K(Tguess,"Nitrogen",pressure_in_words))/(p/101300)
        global K_oxygen = exp(interpolate_K(Tguess,"Oxygen",pressure_in_words))/(p/101300)

        global flash = (frac_N2/(1+(L_on_F*((1/K_nitrogen)-1)))) + ((1-frac_N2)/(1+(L_on_F*((1/K_oxygen)-1))))
        if abs(flash-1) < 0.0005
            break
        elseif flash > 1
            Tguess -= 0.01K
        elseif flash < 1
            Tguess += 0.01K
        end
    end
    Tfinal = Tguess
    x_N2 = (1-K_oxygen)/(K_nitrogen-K_oxygen) #mole fraction of nitrogen in liquid phase
    y_N2 = K_nitrogen*x_N2 #mole fraction of nitrogen in vapor phase

    x_O2 = 1-x_N2
    y_O2 = 1-y_N2
    return Tfinal, y_N2,x_N2,y_O2,x_O2
end