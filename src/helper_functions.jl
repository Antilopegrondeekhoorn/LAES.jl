# Packages
using CoolProp
using Unitful: m, J, kPa, Pa, kg, K, W, °C, l, MW, g, uconvert, °, ustrip, s, kW, hr,mol,kJ
using DataFrames
using CurveFit

########################################################
# State struct
########################################################

mutable struct AirState
    phase::String #gas,liquid or 2phase
    p::Float64
    T::Float64
    h::Float64
    s::Float64
    mdot::Float64
    y_N2::Float64 #gas phase N2 fraction
    x_N2::Float64 #liquid phase N2 fraction
    liquid_fraction::Float64 
    function AirState(phase,p,T,mdot,y_N2,x_N2,liquid_fraction)
        s = new() 
        s.p = p
        s.T = T
        if phase == "gas"
            s.phase = "gas"
            s.y_N2 = y_N2
            y_O2 = 1-y_N2
            s.h = CoolProp.PropsSI("H","T",T,"P|gas",p,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
            s.s = CoolProp.PropsSI("S","T",T,"P|gas",p,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]") 
            s.x_N2 = 0
            s.liquid_fraction = 0
        elseif phase == "liquid"
            s.phase = "liquid"
            s.x_N2 = x_N2
            x_O2 = 1-x_N2
            s.h = CoolProp.PropsSI("H","T",T,"P|liquid",p,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
            s.s = CoolProp.PropsSI("S","T",T,"P|liquid",p,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
            s.y_N2 = 0
            s.liquid_fraction = 1
        elseif phase == "2phase"
            s.phase = "2phase"
            s.y_N2 = y_N2
            y_O2 = 1-y_N2
            s.x_N2 = x_N2
            x_O2 = 1-x_N2
            
            h_liquid = CoolProp.PropsSI("H","P|liquid",p,"T",T,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
            h_vapor = CoolProp.PropsSI("H","P|gas",p,"T",T,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
            s.h = h_liquid*liquid_fraction + h_vapor*(1-liquid_fraction)
    
            s_liquid = CoolProp.PropsSI("S","P|liquid",p,"T",T,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
            s_vapor = CoolProp.PropsSI("S","P|gas",p,"T",T,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
            s.s = s_liquid*liquid_fraction + s_vapor*(1-liquid_fraction)
            s.liquid_fraction = liquid_fraction
        else
            @error("Phase not implemented.")
        end
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

mutable struct CoolantState
    fluid::String
    p::Float64
    T::Float64
    h::Float64
    mdot::Float64
    function CoolantState(fluid,p,T,mdot)
        s = new()
        s.fluid = fluid
        s.p = p
        s.T = T
        s.h = CoolProp.PropsSI("H","T",T,"P",p,fluid)
        s.mdot = mdot
        return s
    end   
end

function State(fluid,p,T,mdot;phase=nothing,y_N2=nothing,x_N2=nothing,liquid_fraction=nothing)
    if fluid == "Air"
        return AirState(phase,p,T,mdot,y_N2,x_N2,liquid_fraction)
    elseif fluid == "Essotherm650"
        return OilState(p,T,mdot)
    elseif fluid == "Methanol" || fluid == "Propane"
        return CoolantState(fluid,p,T,mdot)
    else 
        @error("Fluid not known")
    end
end
    
 
########################################################
# Other helper functions
########################################################

# Compressors and expanders

function isentropic_T(state_in,p_out)
    # a guess of the isentropic temperature at higher/lower pressure
    if ustrip(p_out) > state_in.p
        T_out_is_guess = (state_in.T+10)K
    else
        T_out_is_guess = (state_in.T-10)K
    end
    iterations = 0
    prev_steps = []
    while true
        iterations +=1
        if state_in.phase == "gas"
            s_is_guess = CoolProp.PropsSI("S","T",T_out_is_guess,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.y_N2)]&Oxygen[$(1-state_in.y_N2)]")  #composition doesn't change
        elseif state_in.phase == "liquid"
            s_is_guess = CoolProp.PropsSI("S","T",T_out_is_guess,"P|liquid",ustrip(p_out),"PR::Nitrogen[$(state_in.x_N2)]&Oxygen[$(1-state_in.x_N2)]")  #composition doesn't change
        else
            @error("The calculations of 'isentropic_T' can't be done in the phase of the input state.")
        end
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
        if diff < 0.08
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
                @error("Infinite loop: adjust the steps")
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
        if state_in.phase == "gas"
            h_out_real_guess = CoolProp.PropsSI("H","T",T_out_real_guess,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.y_N2)]&Oxygen[$(1-state_in.y_N2)]")  #composition doesn't change
        elseif state_in.phase == "liquid"
            h_out_real_guess = CoolProp.PropsSI("H","T",T_out_real_guess,"P|liquid",ustrip(p_out),"PR::Nitrogen[$(state_in.x_N2)]&Oxygen[$(1-state_in.x_N2)]")  #composition doesn't change
        else
            @error("The calculations of 'T_real' can't be done in the phase of the input state.")
        end
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
                @error("Infinite loop: adjust the steps")
                break
            end
            prev_steps = []
        end   
    end
    return T_out_real
end

# Intercoolers

function h_Essotherm650(T)
    SpecificHeat = [0.0 1770; 20.0 1850; 40.0 1920; 60.0 1990; 80.0 2060; 100.0 2130; 150.0 2310; 200.0 2490; 250.0 2670; 300.0 2850; 320.0 2920]
    p = poly_fit(SpecificHeat[:,1] + ones(size(SpecificHeat)[1],1)*273.15,SpecificHeat[:,2], 2)
    cp(T) = p[1]+p[2]*T+p[3]*T^2
    h(T) = ((cp(T)+cp(273.15))/2)*(T-273.15)
    return h(T)
end

function T_Essotherm650(h)
    # Rough search
    T = 200:1:1000
    hs = []
    for i = 1:length(T)
        x = h_Essotherm650(T[i])
        push!(hs,x)
    end
    index = findfirst(x->x>h,hs)
    # Fine search
    T = T[index]-1:0.001:T[index]+1
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
        @error("Component of the air not known. Air is considered a mixture of Nitrogen and Oxygen")
    end
    new_value = value1 + (T-T1)/(T2-T1)*(value2-value1)
    return new_value
end

function flash_calculation(p,frac_N2,yield,pressure_in_words) #p given with units
    T_bubble = CoolProp.PropsSI("T","P",p,"Q",0,"PR::Nitrogen[$(frac_N2)]&Oxygen[$(1-frac_N2)]")
    T_dew = CoolProp.PropsSI("T","P",p,"Q",1,"PR::Nitrogen[$(frac_N2)]&Oxygen[$(1-frac_N2)]")

    #Assume 80K (between Tbubble and Tdew)
    T_guess = ((round(ustrip(T_bubble))+round(ustrip(T_dew)))/2)*1.0K
    p = ustrip(p) #[Pa]
    while true
        global L_on_F = yield

        global K_nitrogen = exp(interpolate_K(T_guess,"Nitrogen",pressure_in_words))/(p/101300)
        global K_oxygen = exp(interpolate_K(T_guess,"Oxygen",pressure_in_words))/(p/101300)

        global flash = (frac_N2/(1+(L_on_F*((1/K_nitrogen)-1)))) + ((1-frac_N2)/(1+(L_on_F*((1/K_oxygen)-1))))
        if abs(flash-1) < 0.0005
            break
        elseif flash > 1
            T_guess -= 0.001K
        elseif flash < 1
            T_guess += 0.001K
        end
    end
    T_final = T_guess
    x_N2 = (1-K_oxygen)/(K_nitrogen-K_oxygen)#mole fraction of nitrogen in liquid phase
    y_N2 = K_nitrogen*x_N2 #mole fraction of nitrogen in vapor phase

    x_O2 = 1-x_N2
    y_O2 = 1-y_N2
    return T_final, y_N2,x_N2,y_O2,x_O2
end

# Coldbox
    
function T_compressed_air(h,shifted_hs_compressed_air,T)
    first_diff = shifted_hs_compressed_air[1]-shifted_hs_compressed_air[2]
    if first_diff > 30
        @warn("Make the 'hs_compressed_air' finer. The difference in enthalpy between the first two values should be smaller than 30 J/kg.")
    end
    index = findfirst(x -> x<=h,shifted_hs_compressed_air)
    return T[index]
end

function pinch_coldbox_optimal(state_compressed_air_in,pinch_coldbox,methanol_min,methanol_max,propane_min,propane_max,η_e,pressure_loss_IC)
    # Compressed air
    T_compressed_air_out = 98
    h_airs = []
    T = [state_compressed_air_in.T:-0.02:T_compressed_air_out;T_compressed_air_out]
    for i = 1:length(T)
        compressed_air_h = CoolProp.PropsSI("H","P|gas",state_compressed_air_in.p,"T",T[i],"PR::Nitrogen[$(state_compressed_air_in.y_N2)]&Oxygen[$(1-state_compressed_air_in.y_N2)]")
        push!(h_airs,compressed_air_h)
    end
    h_air_min = minimum(h_airs)
    shifted_h_airs = (h_airs-ones(length(h_airs),1)*h_air_min)
    global h_max_compressed_air = maximum(shifted_h_airs)
    
    # Cold air
    state3_p = state_compressed_air_in.p-state_compressed_air_in.p*pressure_loss_IC
    state4 = State("Air",state3_p-state3_p*pressure_loss_IC,T_compressed_air_out,state_compressed_air_in.mdot;phase = state_compressed_air_in.phase,y_N2 = state_compressed_air_in.y_N2,x_N2 = state_compressed_air_in.x_N2,liquid_fraction = state_compressed_air_in.liquid_fraction)
    state5 =  isentropic_cryoexpander(state4,102000Pa,η_e)
    state6,state7 = separator(state5)

    T_cold_air_in = state5.T
    T_cold_air_out_guess = 284
    iterations = 0
    prev_steps = []
    while true
        iterations +=1
        h_cold_air_in = CoolProp.PropsSI("H","P|gas",102000,"T",T_cold_air_in,"PR::Nitrogen[$(state7.y_N2)]&Oxygen[$(1-state7.y_N2)]")
        h_cold_air_out = CoolProp.PropsSI("H","P|gas",100000,"T",T_cold_air_out_guess,"PR::Nitrogen[$(state7.y_N2)]&Oxygen[$(1-state7.y_N2)]")

        # Coolants
        T_min_propane = propane_min.T
        T_max_propane = propane_max.T
        mdot_propane = propane_min.mdot
        h_min_propane = CoolProp.PropsSI("H","T",T_min_propane,"P|liquid",100000,"Propane")
        h_max_propane = CoolProp.PropsSI("H","T",T_max_propane,"P|liquid",100000,"Propane")
    
        T_min_methanol = methanol_min.T
        T_max_methanol = methanol_max.T
        mdot_methanol = methanol_min.mdot
        h_min_methanol = CoolProp.PropsSI("H","T",T_min_methanol,"P|liquid",100000,"Methanol")
        h_max_methanol = CoolProp.PropsSI("H","T",T_max_methanol,"P|liquid",100000,"Methanol")
    
        #shift propane and methanol
        h_max_propane_shifted = (h_max_propane+abs(h_min_propane))*mdot_propane
        h_max_methanol_shifted = (h_max_methanol+abs(h_min_methanol))*mdot_methanol 
        coeff = h_max_compressed_air/(h_max_propane_shifted+h_max_methanol_shifted)

        h_max_propane_shifted = (h_max_propane+abs(h_min_propane))*mdot_propane*coeff

    
        p = poly_fit(range(h_max_compressed_air,0,1000),range(T_cold_air_out_guess,T_cold_air_in,1000),1)
        Qcoldair(h) = p[1]+p[2]*h

        p2 = poly_fit(range(h_max_propane_shifted,0,1000),range(T_max_propane,T_min_propane,1000),1)
        Qpropane(h) = p2[1]+p2[2]*h

        search_pinch = []

    
        for i = [0:100:h_max_compressed_air;h_max_compressed_air]
            difference = T_compressed_air(i,shifted_h_airs,T)-Qcoldair(i)
            push!(search_pinch,difference)
        end
        pinch = minimum(search_pinch)
        
        diff = abs(pinch_coldbox - pinch)
        if diff > 5
            step = 4
        elseif 2 < diff < 5
            step = 2
        elseif 1 < diff < 2
            step = 1
        elseif 0.5 < diff < 1
            step = 0.5
        elseif 0.2 < diff < 0.5
            step = 0.2
        elseif 0.1 < diff < 0.2
            step = 0.1
        elseif 0.02 < diff < 0.1
            step = 0.05
        elseif 0.005 < diff < 0.02
            step = 0.005
        else
            step = 0.001
        end

        if diff < 0.001
            global T_cold_air_out = T_cold_air_out_guess
            break
        elseif pinch < pinch_coldbox
            T_cold_air_out_guess -= step
            push!(prev_steps,-1)
        elseif pinch > pinch_coldbox
            T_cold_air_out_guess += step
            push!(prev_steps,+1)
        end

        #prevent being in an infinite loop
        if iterations%4 == 0
            if sum(prev_steps) ==  0
                @error("Infinite loop: adjust the steps")
                break
            end
            prev_steps = []
        end 
    end
    return state6.mdot,T_cold_air_out
end

function pinch_coldbox_p_less_optimal(state_compressed_air_in,pinch_coldbox,T_cold_air_out_optimal,methanol_min,methanol_max,propane_min,propane_max,η_e,pressure_loss)
    iterations = 0
    prev_steps = []
    T_guess_compressed_air_out = 115
    while true
        iterations +=1
        # Compressed air
        h_airs = []
        T = [state_compressed_air_in.T:-0.02:T_guess_compressed_air_out;T_guess_compressed_air_out]
        for i = 1:length(T)
            compressed_air_h = CoolProp.PropsSI("H","P|gas",state_compressed_air_in.p,"T",T[i],"PR::Nitrogen[$(state_compressed_air_in.y_N2)]&Oxygen[$(1-state_compressed_air_in.y_N2)]")
            push!(h_airs,compressed_air_h)
        end
        h_air_min = minimum(h_airs)
        shifted_h_airs = (h_airs-ones(length(h_airs),1)*h_air_min)
        global h_max_compressed_air = maximum(shifted_h_airs)
    
        # Cold air
        state3_p = state_compressed_air_in.p-state_compressed_air_in.p*pressure_loss
        state4 = State("Air",state3_p-state3_p*pressure_loss,T_guess_compressed_air_out,state_compressed_air_in.mdot;phase = state_compressed_air_in.phase,y_N2 = state_compressed_air_in.y_N2,x_N2 = state_compressed_air_in.x_N2,liquid_fraction = state_compressed_air_in.liquid_fraction)
        state5 =  isentropic_cryoexpander(state4,102000Pa,η_e)
        global state6,state7 = separator(state5)

        T_cold_air_in = state5.T 
        h_cold_air_in = CoolProp.PropsSI("H","P|gas",102000,"T",T_cold_air_in,"PR::Nitrogen[$(state7.y_N2)]&Oxygen[$(1-state7.y_N2)]")
        h_cold_air_out = CoolProp.PropsSI("H","P|gas",100000,"T",T_cold_air_out_optimal,"PR::Nitrogen[$(state7.y_N2)]&Oxygen[$(1-state7.y_N2)]") #from the optimal case

        p = poly_fit(range(h_max_compressed_air,0,1000),range(T_cold_air_out_optimal,T_cold_air_in,1000),1)
        Qcoldair(h) = p[1]+p[2]*h

        search_pinch = []
    
        for i = [0:100:h_max_compressed_air;h_max_compressed_air]
            difference = T_compressed_air(i,shifted_h_airs,T)-Qcoldair(i)
            push!(search_pinch,difference)
        end
        pinch = minimum(search_pinch)
    
        diff = abs(pinch_coldbox - pinch)
        
        if diff > 5
            step = 5
        elseif 2 < diff < 5
            step = 2
        elseif 0.5 < diff < 2
            step = 1
        elseif 0.2 < diff < 0.5
            step = 0.5
        elseif 0.1 < diff < 0.2
            step = 0.1
        elseif 0.02 < diff < 0.1
            step = 0.05
        elseif 0.01 < diff < 0.02
            step = 0.002
        else
            step = 0.001
        end
        
        if diff < 0.005
            global T_compressed_air_out = T_guess_compressed_air_out
            break
        elseif pinch < pinch_coldbox
            T_guess_compressed_air_out += step
            push!(prev_steps,+1)
        elseif pinch > pinch_coldbox
            T_guess_compressed_air_out -= step
            push!(prev_steps,-1)
        end
        
        #prevent being in an infinite loop
        if iterations%4 == 0
            if sum(prev_steps) ==  0
                @error("Infinite loop: adjust the steps")
                break
            end
            prev_steps = []
        end 
    end

    return state6.mdot,T_compressed_air_out
end

function pinch_coldbox_p_more_optimal(state_compressed_air_in,pinch_coldbox,T_compressed_air_out_optimal,methanol_min,methanol_max,propane_min,propane_max,η_e,pressure_loss)
    iterations = 0
    prev_steps = []
    T_guess_cold_air_out = 300
    while true
        iterations +=1
        # Compressed air
        h_airs = []
        T = [state_compressed_air_in.T:-0.02:T_compressed_air_out_optimal;T_compressed_air_out_optimal]
        for i = 1:length(T)
            compressed_air_h = CoolProp.PropsSI("H","P|gas",state_compressed_air_in.p,"T",T[i],"PR::Nitrogen[$(state_compressed_air_in.y_N2)]&Oxygen[$(1-state_compressed_air_in.y_N2)]")
            push!(h_airs,compressed_air_h)
        end
        h_air_min = minimum(h_airs)
        shifted_h_airs = (h_airs-ones(length(h_airs),1)*h_air_min)
        global h_max_compressed_air = maximum(shifted_h_airs)
    
        # Cold air
        state3_p = state_compressed_air_in.p-state_compressed_air_in.p*pressure_loss
        state4 = State("Air",state3_p-state3_p*pressure_loss,T_compressed_air_out_optimal,state_compressed_air_in.mdot;phase = state_compressed_air_in.phase,y_N2 = state_compressed_air_in.y_N2,x_N2 = state_compressed_air_in.x_N2,liquid_fraction = state_compressed_air_in.liquid_fraction)
        state5 =  isentropic_cryoexpander(state4,102000Pa,η_e)
        global state6,state7 = separator(state5)

        T_cold_air_in = state5.T 
        h_cold_air_in = CoolProp.PropsSI("H","P|gas",102000,"T",T_cold_air_in,"PR::Nitrogen[$(state7.y_N2)]&Oxygen[$(1-state7.y_N2)]")
        h_cold_air_out = ((state_compressed_air_in.mdot*(state_compressed_air_in.h-state4.h)- propane_min.mdot*(propane_max.h-propane_min.h) - methanol_min.mdot*(methanol_max.h-methanol_min.h))/state7.mdot)+state7.h

        h_cold_air_out_guess = CoolProp.PropsSI("H","P|gas",100000,"T",T_guess_cold_air_out,"PR::Nitrogen[$(state7.y_N2)]&Oxygen[$(1-state7.y_N2)]") 
        diff = abs(h_cold_air_out- h_cold_air_out_guess)
        
        if diff > 10000
            step = 10
        elseif 5000 < diff < 10000
            step = 5
        elseif 2000 < diff < 5000
            step = 2
        elseif 1000 < diff < 2000
            step = 1
        elseif 500 < diff < 1000
            step = 0.5
        elseif 200 < diff < 500
            step = 0.1
        elseif 50 < diff < 200
            step = 0.05
        elseif 20 < diff < 50
            step = 0.01
        elseif 10 < diff < 20
            step = 0.005
        elseif 5 < diff < 10
            step = 0.002
        else
            step = 0.001
        end
        
        if diff < 5
            global T_cold_air_out = T_guess_cold_air_out
            break
        elseif h_cold_air_out_guess < h_cold_air_out
            T_guess_cold_air_out += step
            push!(prev_steps,+1)
        elseif h_cold_air_out_guess  > h_cold_air_out
            T_guess_cold_air_out -= step
            push!(prev_steps,-1)
        end

        #prevent being in an infinite loop
        if iterations%4 == 0
            if sum(prev_steps) ==  0
                @error("Infinite loop: adjust the steps")
                break
            end
            prev_steps = []
        end 
    end

    return state6.mdot,T_cold_air_out
end

# Extra function for the charging cycle
function compare(sol1::AirState,sol2::AirState)
    if sol1.p == sol2.p && sol1.T == sol2.T && sol1.mdot == sol2.mdot && sol1.y_N2 == sol2.y_N2 && sol1.x_N2 == sol2.x_N2 && sol1.phase == sol2.phase
        return true
    else
        return false
    end   
end

# Extra function for the discharge cycle

function find_T_5R(h_5R,p_5R,y_N2)
    T_5R_guess = 430K  
    iterations = 0
    prev_steps = []
    while true
        iterations +=1
        h_5R_guess = CoolProp.PropsSI("H","T",T_5R_guess,"P|gas",p_5R,"PR::Nitrogen[$(y_N2)]&Oxygen[$(1-y_N2)]") 
        #change step when closer to solution
        diff = abs(ustrip(h_5R_guess)-h_5R) 
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
            global T_5R = T_5R_guess
            break
        elseif ustrip(h_5R_guess) < h_5R
            T_5R_guess += step
            push!(prev_steps,+1)
        elseif ustrip(h_5R_guess) > h_5R
            T_5R_guess -= step
            push!(prev_steps,-1)
        end
        #prevent being in an infinite loop
        if iterations%4 == 0
            if sum(prev_steps) ==  0
                @error("Infinite loop: adjust the steps")
                break
            end
            prev_steps = []
        end   
    end 
    return T_5R
end
     