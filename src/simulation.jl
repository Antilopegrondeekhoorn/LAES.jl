include("components.jl")

function charging_cycle(state_in,pinch_IC,pinch_coldbox,pressure_loss_IC,stateOil_in,oil_distribution,methanol_min,methanol_max,propane_min,propane_max,η_c,η_e)
    solutions = []
    state1 = state_in
    #distribution calculated from reference paper
    stateOil_in1 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[1])
    stateOil_in2 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[2])
    
    for i = 1:20
        global state2A = isentropic_compressor(state1,1480000Pa,η_c)
        global state2B,stateOil_out1 = intercooler("Cool",state2A,stateOil_in1,pinch_IC,pressure_loss_IC)
        global state2C = isentropic_compressor(state2B,18098000Pa,η_c)
        global state2,stateOil_out2 = intercooler("Cool",state2C,stateOil_in2,pinch_IC,pressure_loss_IC)
        state3_p = state2.p-state2.p*pressure_loss_IC
        global state4 = State("Air",state3_p-state3_p*pressure_loss_IC,propane_min.T+pinch_coldbox,state2.mdot;phase = state2.phase,y_N2 = state2.y_N2,x_N2 = state2.x_N2,liquid_fraction = state2.liquid_fraction)
    
        global state5 =  isentropic_cryoexpander(state4,102000Pa,η_e)
        global state6,state7 = separator(state5)
        global state9 = coldbox(methanol_min,methanol_max,propane_min,propane_max,state2,state4,state7)
        state9.p = p10
        global state10 = State("Air",p10,T10,state6.mdot;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = state9.liquid_fraction)#standard conditions
        global state1 = State("Air",p10,state9.mdot*state9.T+state10.mdot*state10.T,mdot1;phase = "gas",y_N2 = state9.mdot*state9.y_N2+state10.mdot*state10.y_N2,x_N2 = x_N2,liquid_fraction = state10.liquid_fraction)
        #println(i)
        #println(state1) #check convergence
        push!(solutions,state1)
        if length(solutions) > 1 && compare(solutions[end-1],solutions[end])
            break
        end
        if i == 20
            @warn("Did not converge after 20 iterations")
        end
    end
    return state1,state2A,state2B,state2C,state2,state4,state5,state6,state7,state9,state10,stateOil_out1,stateOil_out2
end

function discharge(state1R,p_2R,stateOil_in,oil_distribution,propane_max,methanol_max,pinch_coldbox,pinch_superheaters,η_e,η_pump,pressure_loss)
    state2R = isentropic_cryopump(state1R,p_2R,η_pump)
    state3R = heater_coldstorage(state2R,propane_max,pinch_coldbox,pressure_loss)
    state4R = heater_coldstorage(state3R,methanol_max,pinch_coldbox,pressure_loss)
    
    #distribution calculated from reference paper
    stateOil_in1 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[1])
    stateOil_in2 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[2])
    stateOil_in3 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[3])
        
    #calculate discharge process with guessed state 5R
    state_5R_guess = State("Air",p_5R,T_guess_5R,state4R.mdot;phase = "gas",y_N2 = state4R.y_N2,x_N2 = state4R.x_N2,liquid_fraction = state4R.liquid_fraction)
    #discharging
    state6R,stateOil_out1 = intercooler("Heat",state_5R_guess,stateOil_in1,pinch_superheaters,pressure_loss)
    state7R = isentropic_expander(state6R,1590000,η_e)
    state8R,stateOil_out2 = intercooler("Heat",state7R,stateOil_in2,pinch_superheaters,pressure_loss)
    state9R = isentropic_expander(state8R,401000,η_e)
    state10R,stateOil_out3 = intercooler("Heat",state9R,stateOil_in3,pinch_superheaters,pressure_loss)
    state11R = isentropic_expander(state10R,101000,η_e)
    state12R = State("Air",state11R.p-state11R.p*pressure_loss,state4R.T+pinch_coldbox,state11R.mdot;phase = state11R.phase,y_N2 = state11R.y_N2,x_N2 = state11R.x_N2,liquid_fraction = state11R.liquid_fraction)
    
    #energy balance over heat exchanger
    h_5R = state11R.h-state12R.h+state4R.h
    T_5R = find_T_5R(h_5R,state_5R_guess.y_N2)
        
    #define state5R
    state5R = State("Air",p_5R,ustrip(T_5R),state4R.mdot;phase = "gas",y_N2 = state4R.y_N2,x_N2 = state4R.x_N2,liquid_fraction = state4R.liquid_fraction)
    #discharging
    global state6R,stateOil_out1 = intercooler("Heat",state_5R_guess,stateOil_in1,pinch_superheaters,pressure_loss)
    global state7R = isentropic_expander(state6R,1590000,η_e)
    global state8R,stateOil_out2 = intercooler("Heat",state7R,stateOil_in2,pinch_superheaters,pressure_loss)
    global state9R = isentropic_expander(state8R,401000,η_e)
    global state10R,stateOil_out3 = intercooler("Heat",state9R,stateOil_in3,pinch_superheaters,pressure_loss)
    global state11R = isentropic_expander(state10R,101000,η_e)
    global state12R = State("Air",state11R.p-state11R.p*pressure_loss,state4R.T+pinch_coldbox,state11R.mdot;phase = state11R.phase,y_N2 = state11R.y_N2,x_N2 = state11R.x_N2,liquid_fraction = state11R.liquid_fraction)
    return state1R,state2R,state3R,state4R,state5R,state6R,state7R,state8R,state9R,state10R,state11R,state12R,stateOil_out1,stateOil_out2,stateOil_out3
end