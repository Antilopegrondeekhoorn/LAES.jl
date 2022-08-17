include("helper_functions.jl")

#########################################################################
# Charging cycle
#########################################################################

function isentropic_compressor(state_in::AirState,p_out,η_c)
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
    iterations = 0
    prev_steps = []
    yield_guess_is = 0.88    
    while true
        iterations +=1
        T_out, y_N2,x_N2,y_O2,x_O2 = flash_calculation(p_out,state_in.y_N2,yield_guess_is,"1atm") 

        h_liquid = CoolProp.PropsSI("H","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        h_vapor = CoolProp.PropsSI("H","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        global h_out_is = h_liquid*yield_guess_is + h_vapor*(1-yield_guess_is)
        
        s_liquid = CoolProp.PropsSI("S","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        s_vapor = CoolProp.PropsSI("S","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        s_out_is = s_liquid*yield_guess_is + s_vapor*(1-yield_guess_is)
        
        diff = abs(state_in.s-ustrip(s_out_is))
       
        if diff > 300
            step = 0.1
        elseif 150 < diff < 300
            step = 0.01
        elseif 75 < diff < 150
            step = 0.005
        elseif 25 < diff < 75
            step = 0.001
        elseif 10 < diff < 25
            step = 0.0005
        elseif 4 < diff < 10
            step = 0.0001
        else
            step = 0.00005
        end
        
        if diff <0.5 #arbitrary value
            global yield_is = yield_guess_is
            break
        elseif ustrip(s_out_is) < state_in.s
            yield_guess_is -=step
            push!(prev_steps,-1)
        elseif ustrip(s_out_is) > state_in.s
            yield_guess_is +=step
            push!(prev_steps,+1)
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
    h_out_real = state_in.h- η_e*(state_in.h-ustrip(h_out_is))
    
    # Real calculations
    iterations = 0
    prev_steps = []
    yield_guess_real = 0.84
    while true
        iterations +=1
        global T_out, y_N2,x_N2,y_O2,x_O2 = flash_calculation(p_out,state_in.y_N2,yield_guess_real,"1atm") 
        
        h_liquid = CoolProp.PropsSI("H","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        h_vapor = CoolProp.PropsSI("H","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        global h_out = h_liquid*yield_guess_real + h_vapor*(1-yield_guess_real)
        
        s_liquid = CoolProp.PropsSI("S","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        s_vapor = CoolProp.PropsSI("S","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        global s_out = s_liquid*yield_guess_real + s_vapor*(1-yield_guess_real)
        
        #change step when closer to solution
        diff = abs(h_out_real-ustrip(h_out))
        
        if diff > 5000
            step = 0.02
        elseif 2500 < diff < 5000
            step = 0.01
        elseif 1000 < diff < 2500
            step = 0.005
        elseif 400 < diff < 1000
            step = 0.001
        elseif 25 < diff < 400
            step = 0.0001
        else
            step = 0.00005
        end
        
        #change the yield guess
        if diff < 15
            global yield_real = yield_guess_real
            break
        elseif ustrip(h_out) < h_out_real
            yield_guess_real -=step
            push!(prev_steps,-1)
        elseif ustrip(h_out) > h_out_real
            yield_guess_real +=step
            push!(prev_steps,+1)
        end
        
        #prevent being in an infinite loop
        if iterations%8 == 0
            if sum(prev_steps) ==  0
                @warn "Infinite loop: adjust the steps"
                break
            end
            prev_steps = []
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

#########################################################################
# Storage
#########################################################################

function storage_tank(state_in,p_storage)
    T1R = CoolProp.PropsSI("T","P",p_storage,"Q",0,"PR::Nitrogen[$(state_in.x_N2)]&Oxygen[$(1-state_in.x_N2)]") #stored as saturated liquid
    state_out = State("Air",100000,ustrip(T1R),state_in.mdot;phase = "liquid",y_N2 = state_in.y_N2,x_N2 = state_in.x_N2,liquid_fraction = state_in.liquid_fraction)
    return state_out
end

#########################################################################
# Discharge cycle
#########################################################################

function isentropic_cryopump(state_in::AirState,p_out,η_pump)
    T_out_is = isentropic_T(state_in,p_out)
    h_out_is = CoolProp.PropsSI("H","T",T_out_is,"P|liquid",ustrip(p_out),"PR::Nitrogen[$(state_in.x_N2)]&Oxygen[$(1-state_in.x_N2)]")
    h_out_real = (ustrip(h_out_is)-state_in.h)/η_pump + state_in.h
    T_out_real = T_real(T_out_is,state_in,p_out,h_out_real)
    
    state_out = State("Air",ustrip(p_out),ustrip(T_out_real),state_in.mdot;phase = "gas",y_N2 = state_in.x_N2,x_N2 =0.0,liquid_fraction = 0.0)
    return state_out
end

function heater_coldstorage(state_in::AirState,coolant_max::CoolantState,pinch,pressure_loss)
    state_out = State("Air",state_in.p-state_in.p*pressure_loss,coolant_max.T-pinch,state_in.mdot;phase = state_in.phase,y_N2 = state_in.y_N2,x_N2 = state_in.x_N2,liquid_fraction = state_in.liquid_fraction)
    return state_out
end

function isentropic_expander(state_in::AirState,p_out,η_e)
    T_out_is = isentropic_T(state_in,p_out)
    h_out_is = CoolProp.PropsSI("H","T",T_out_is,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.y_N2)]&Oxygen[$(1-state_in.y_N2)]")
    h_out_real = state_in.h- η_e*(state_in.h-ustrip(h_out_is))
    T_out_real = T_real(T_out_is,state_in,p_out,h_out_real)
    state_out = State("Air",ustrip(p_out),ustrip(T_out_real),state_in.mdot;phase = "gas",y_N2 = state_in.y_N2,x_N2 = state_in.x_N2,liquid_fraction = state_in.liquid_fraction)
    return state_out
end   