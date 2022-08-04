include("helper_functions.jl")

function isentropic_compressor(state_in::AirState,p_out,η_c)
    T_out_is = isentropic_T(state_in,p_out)
    h_out_is = CoolProp.PropsSI("H","T",T_out_is,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.frac_N2)]&Oxygen[$(1-state_in.frac_N2)]")
    h_out_real = (ustrip(h_out_is)-state_in.h)/η_c + state_in.h
    T_out_real = T_real(T_out_is,state_in,p_out,h_out_real)
    
    state_out = State("Air",ustrip(p_out),ustrip(T_out_real),state_in.mdot;frac_N2 = state_in.frac_N2)
    return state_out
end

function intercooler(state_in_air::AirState,state_in_oil::OilState,pinch_IC,pressureloss_IC)
    #pinch at the cold end of the intercooler???
    T_out_air = state_in_oil.T+pinch_IC
    #pressure loss
    p_out_air = state_in_air.p - state_in_air.p*pressureloss_IC
    #output state
    state_out_air = State("Air",p_out_air,T_out_air,state_in_air.mdot;frac_N2 = state_in_air.frac_N2)
    
    #energy balance
    h_out_oil = (2*(state_in_air.mdot/state_in_oil.mdot)*(state_in_air.h-state_out_air.h)) + state_in_oil.h
    T_out_oil = T_Essotherm650(h_out_oil)
    state_out_oil = State("Essotherm650",100000,T_out_oil,state_in_oil.mdot)
    
    return state_out_air,state_out_oil
end

function separator(state_in,yield)
    Tout, y_N2,x_N2,y_O2,x_O2 = flash_calculation(state_in.p,state_in.frac_N2,yield,"1atm") 

    h_liquid = CoolProp.PropsSI("H","P|liquid",state_in.p,"T",Tout,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
    h_vapor = CoolProp.PropsSI("H","P|gas",state_in.p,"T",Tout,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
    hout = h_liquid*yield + h_vapor*(1-yield)
    
    s_liquid = CoolProp.PropsSI("S","P|liquid",state_in.p,"T",Tout,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
    s_vapor = CoolProp.PropsSI("S","P|gas",state_in.p,"T",Tout,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
    sout = s_liquid*yield + s_vapor*(1-yield)
    
    liquidstate = State("Air",state_in.p,ustrip(Tout),state2.mdot*yield;frac_N2 = x_N2) #wrong calculation of h and s
    liquidstate.h = ustrip(h_liquid)
    liquidstate.s = ustrip(s_liquid)

    vaporstate = State("Air",state_in.p,ustrip(Tout),state2.mdot*(1-yield);frac_N2 = y_N2) #wrong calculation of h and s
    vaporstate.h = ustrip(h_vapor)
    vaporstate.s = ustrip(s_vapor)
    
    return liquidstate,vaporstate
end

function isentropic_cryoexpander(state_in::AirState,p_out,η_e)
    #isentropic calculations
    yield_guess_is = 0.88
    while true
        T_out, y_N2,x_N2,y_O2,x_O2 = flash_calculation(p_out,state_in.frac_N2,yield_guess_is,"1atm") 

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
        global T_out, y_N2,x_N2,y_O2,x_O2 = flash_calculation(p_out,state_in.frac_N2,yield_guess_real,"1atm") 
        
        h_liquid = CoolProp.PropsSI("H","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        h_vapor = CoolProp.PropsSI("H","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        global h_out = h_liquid*yield_guess_real + h_vapor*(1-yield_guess_real)
        
        s_liquid = CoolProp.PropsSI("S","P|liquid",ustrip(p_out),"T",T_out,"PR::Nitrogen[$x_N2]&Oxygen[$x_O2]")
        s_vapor = CoolProp.PropsSI("S","P|gas",ustrip(p_out),"T",T_out,"PR::Nitrogen[$y_N2]&Oxygen[$y_O2]")
        global s_out = s_liquid*yield_guess_real + s_vapor*(1-yield_guess_real)

        if abs(h_out_real-ustrip(h_out)) <25
            global yield = yield_guess_real
            break
        elseif ustrip(h_out) < h_out_real
            yield_guess_real -=0.0001
        elseif ustrip(h_out) > h_out_real
            yield_guess_real +=0.0001
        end
    end
    state_out = State("Air",ustrip(p_out),ustrip(T_out),state_in.mdot;frac_N2=y_N2)
    state_out.h = ustrip(h_out) #State doesn't calculate in 2phase domain
    state_out.s = ustrip(s_out)
    return state_out,yield
end