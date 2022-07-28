include("helper_functions.jl")

function isentropic_compressor(state_in,p_out,η_c)
    T_out_is = isentropic_T(state_in,p_out)
    h_out_is = CoolProp.PropsSI("H","T",T_out_is,"P|gas",ustrip(p_out),"PR::Nitrogen[$(state_in.frac_N2)]&Oxygen[$(1-state_in.frac_N2)]")
    h_out_real = (ustrip(h_out_is)-state_in.h)/η_c + state_in.h
    T_out_real = T_real(T_out_is,state_in,p_out,h_out_real)
    
    state_out = State(ustrip(p_out),ustrip(T_out_real),state_in.frac_N2,state_in.mdot)
    return state_out
end

