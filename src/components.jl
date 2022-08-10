include("helper_functions.jl")

#########################################################################
# Charging cycle
#########################################################################

function isentropic_compressor(state_in,p_out,η_c)
    T_out_is = isentropic_T(state_in,p_out)
    h_out_is = CoolProp.PropsSI("H","T",T_out_is,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.y_N2)]&Oxygen[$(1-state_in.y_N2)]")
    h_out_real = (ustrip(h_out_is)-state_in.h)/η_c + state_in.h
    T_out_real = T_real(T_out_is,state_in,p_out,h_out_real)
    
    state_out = State("Air",ustrip(p_out),ustrip(T_out_real),state_in.mdot;phase = state_in.phase,y_N2 = state_in.y_N2,x_N2 = state_in.x_N2,liquid_fraction = state_in.liquid_fraction)
    return state_out
end

function intercooler(func::String,state_in_air::AirState,state_in_oil::OilState,pinch,pressureloss)
    #pinch at the cold end of the intercooler???
    if func == "Cool"
        T_out_air = state_in_oil.T+pinch
    elseif func == "Heat"
        T_out_air = state_in_oil.T-pinch
    else
        @warn("Function of the intercooler not defined. The intercooler can 'Cool' or 'Heat' the air stream.")
    end
    #pressure loss
    p_out_air = state_in_air.p - state_in_air.p*pressureloss
    #output state
    state_out_air = State("Air",p_out_air,T_out_air,state_in_air.mdot;phase = state_in_air.phase,y_N2 = state_in_air.y_N2,x_N2 = state_in_air.x_N2,liquid_fraction = state_in_air.liquid_fraction)
    
    #energy balance 
    h_out_oil = ((state_in_air.mdot/state_in_oil.mdot)*(state_in_air.h-state_out_air.h)) + state_in_oil.h
    T_out_oil = T_Essotherm650(h_out_oil)
    state_out_oil = State("Essotherm650",100000,T_out_oil,state_in_oil.mdot)
    
    return state_out_air,state_out_oil
end


function isentropic_cryoexpander(state_in::AirState,p_out,η_e)
    #isentropic calculations
    yield_guess_is = 0.88
    while true
        T_out, y_N2,x_N2,y_O2,x_O2 = flash_calculation(p_out,state_in.y_N2,yield_guess_is,"1atm") 

        h_liquid = CoolProp.PropsSI("H","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        h_vapor = CoolProp.PropsSI("H","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        global h_out_is = h_liquid*yield_guess_is + h_vapor*(1-yield_guess_is)
        
        s_liquid = CoolProp.PropsSI("S","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        s_vapor = CoolProp.PropsSI("S","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        s_out_is = s_liquid*yield_guess_is + s_vapor*(1-yield_guess_is)
        
        if abs(state_in.s-ustrip(s_out_is)) <0.2 #arbitrary value
            global yield_is = yield_guess_is
            break
        elseif ustrip(s_out_is) < state_in.s
            yield_guess_is -=0.00005
        elseif ustrip(s_out_is) > state_in.s
            yield_guess_is +=0.00005
        end
    end
    h_out_real = state_in.h- η_e*(state_in.h-ustrip(h_out_is))
    
    # Real calculations
    yield_guess_real = 0.854
    while true
        global T_out, y_N2,x_N2,y_O2,x_O2 = flash_calculation(p_out,state_in.y_N2,yield_guess_real,"1atm") 
        
        h_liquid = CoolProp.PropsSI("H","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        h_vapor = CoolProp.PropsSI("H","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        global h_out = h_liquid*yield_guess_real + h_vapor*(1-yield_guess_real)
        
        s_liquid = CoolProp.PropsSI("S","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        s_vapor = CoolProp.PropsSI("S","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        global s_out = s_liquid*yield_guess_real + s_vapor*(1-yield_guess_real)

        if abs(h_out_real-ustrip(h_out)) <25
            global yield_real = yield_guess_real
            break
        elseif ustrip(h_out) < h_out_real
            yield_guess_real -=0.0001
        elseif ustrip(h_out) > h_out_real
            yield_guess_real +=0.0001
        end
    end
    state_out = State("Air",ustrip(p_out),ustrip(T_out),state_in.mdot;phase = "2phase",y_N2 = y_N2,x_N2 = x_N2,liquid_fraction = yield_real)
    return state_out
end

function separator(state_in::AirState)
    liquidstate = State("Air",state_in.p,state_in.T,state_in.liquid_fraction;phase = "liquid",y_N2 = 0,x_N2 = state_in.x_N2,liquid_fraction = 1) #wrong calculation of h and s
    vaporstate = State("Air",state_in.p,state_in.T,1-state_in.liquid_fraction;phase = "gas",y_N2 = state_in.y_N2,x_N2 = 0,liquid_fraction = 0) #wrong calculation of h and s
    return liquidstate,vaporstate
end

function coldbox(coolant1_min::CoolantState,coolant1_max::CoolantState,coolant2_min::CoolantState,coolant2_max::CoolantState,state_in_air::AirState,state_out_air::AirState,state_in_air_feedback::AirState) #to calculate statecoldout
    #energy balance
    h_air_out_feedback = ((state_in_air.mdot*(state_in_air.h-state_out_air.h) - coolant1_min.mdot*(coolant1_max.h-coolant1_min.h) - coolant2_min.mdot*(coolant2_max.h-coolant2_min.h))/state_in_air_feedback.mdot)+state_in_air_feedback.h
    p_out_intermediate = state_in_air_feedback.p-state_in_air_feedback.p*0.01 #2 times 1% procent, because 2 HEX
    p_out = p_out_intermediate-p_out_intermediate*0.01
    
    T_air_out_feedback_guess = 300.001K #[K]
    iterations = 0
    prev_steps = []
    while true
        iterations +=1
        h_guess = CoolProp.PropsSI("H","P|gas",p_out,"T",T_air_out_feedback_guess,"PR::Nitrogen[$(state_in_air_feedback.y_N2)]&Oxygen[$(1-state_in_air_feedback.y_N2)]")
        
        #change step when closer to solution
        diff = abs(h_air_out_feedback-ustrip(h_guess))
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
            global T_air_out_feedback = T_air_out_feedback_guess
            break
        elseif ustrip(h_guess) < h_air_out_feedback
            T_air_out_feedback_guess += step
            push!(prev_steps,+1)
        elseif ustrip(h_guess) > h_air_out_feedback
            T_air_out_feedback_guess -= step
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
    state_out = State("Air",p_out,ustrip(T_air_out_feedback),state_in_air_feedback.mdot;phase = state_in_air_feedback.phase,y_N2 = state_in_air_feedback.y_N2,x_N2 = state_in_air_feedback.x_N2,liquid_fraction = state_in_air_feedback.liquid_fraction)
    return state_out
end

#########################################################################
# Discharge cycle
#########################################################################

function isentropic_cryopump(state_in,p_out,η_pump)
    T_out_is = isentropic_T(state_in,p_out)
    h_out_is = CoolProp.PropsSI("H","T",T_out_is,"P|liquid",ustrip(p_out),"PR::Nitrogen[$(state_in.x_N2)]&Oxygen[$(1-state_in.x_N2)]")
    h_out_real = (ustrip(h_out_is)-state_in.h)/η_pump + state_in.h
    T_out_real = T_real(T_out_is,state_in,p_out,h_out_real)
    
    state_out = State("Air",ustrip(p_out),ustrip(T_out_real),state_in.mdot;phase = "gas",y_N2 = state_in.x_N2,x_N2 =0.0,liquid_fraction = 0.0)
    return state_out
end

function heater_coldstorage(state_in::AirState,coolant_max::CoolantState,pinch,pressure_loss)
    state_out = State("Air",state_in.p+state_in.p*pressure_loss,coolant_max.T-pinch,state_in.mdot;phase = state_in.phase,y_N2 = state_in.y_N2,x_N2 = state_in.x_N2,liquid_fraction = state_in.liquid_fraction)
    return state_out
end

function isentropic_expander(state_in,p_out,η_e)
    T_out_is = isentropic_T(state_in,p_out)
    h_out_is = CoolProp.PropsSI("H","T",T_out_is,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.y_N2)]&Oxygen[$(1-state_in.y_N2)]")
    h_out_real = state_in.h- η_e*(state_in.h-ustrip(h_out_is))
    T_out_real = T_real(T_out_is,state_in,p_out,h_out_real)
    state_out = State("Air",ustrip(p_out),ustrip(T_out_real),state_in.mdot;phase = "gas",y_N2 = state_in.y_N2,x_N2 = state_in.x_N2,liquid_fraction = state_in.liquid_fraction)
    return state_out
end   