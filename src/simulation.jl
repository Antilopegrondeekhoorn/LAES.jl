include("components.jl")

function charging_cycle(state_in::AirState,ambient_state::AirState,stateOil_in::OilState,oil_distribution,methanol_min,methanol_max,propane_min,propane_max,pinch_IC,pinch_coldbox,compressor_pressures,η_c,η_e,pressure_loss)
    number_of_compressors = 2
    if length(compressor_pressures) != number_of_compressors && length(oil_distribution) != number_of_compressors
        error("The length of the 'compressor_pressures' and the 'oil_distribution' must be equal to the number of compressors.")
    end
    
    solutions = []
    global state1 = state_in

    #distribution calculated from reference paper
    stateOil_in1 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[1])
    stateOil_in2 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[2])
    
    for i = 1:20
        global state2A = isentropic_compressor(state1,compressor_pressures[1],η_c)
        global state2B,stateOil_out1 = intercooler("Cool",state2A,stateOil_in1,pinch_IC,pressure_loss)
        global state2C = isentropic_compressor(state2B,compressor_pressures[2],η_c)
        global state2,stateOil_out2 = intercooler("Cool",state2C,stateOil_in2,pinch_IC,pressure_loss)
        state3_p = state2.p-state2.p*pressure_loss

        if compressor_pressures[end] == 17917000
            global T4 = propane_min.T+pinch_coldbox
            global yield,T9 = pinch_coldbox_optimal(state2,pinch_coldbox,methanol_min,methanol_max,propane_min,propane_max,η_e,pressure_loss)
        elseif compressor_pressures[end] < 17917000 #Assume that the optimal pressure in the reference paper is the effective optimal pressure
            global ~,T9 = pinch_coldbox_optimal(state2,pinch_coldbox,methanol_min,methanol_max,propane_min,propane_max,η_e,pressure_loss)
            global yield,T4 = pinch_coldbox_p_less_optimal(state2,pinch_coldbox,T9,methanol_min,methanol_max,propane_min,propane_max,η_e,pressure_loss)
        elseif compressor_pressures[end] > 17917000
            global T4 = propane_min.T+pinch_coldbox
            global yield,T9 = pinch_coldbox_p_more_optimal(state2,pinch_coldbox,T4,methanol_min,methanol_max,propane_min,propane_max,η_e,pressure_loss)
        end

        global state4 = State("Air",state3_p-state3_p*pressure_loss,T4,state2.mdot;phase = state2.phase,y_N2 = state2.y_N2,x_N2 = state2.x_N2,liquid_fraction = state2.liquid_fraction)
        global state5 =  isentropic_cryoexpander(state4,102000,η_e)
        global state6,state7 = separator(state5)
        
        global state9 = State("Air",state1.p,T9,state7.mdot;phase = state7.phase,y_N2 = state7.y_N2,x_N2 = state7.x_N2,liquid_fraction = state7.liquid_fraction)
        global state10 = State("Air",ambient_state.p,ambient_state.T,state6.mdot;phase = ambient_state.phase,y_N2 = ambient_state.y_N2,x_N2 = ambient_state.x_N2,liquid_fraction = ambient_state.liquid_fraction) #standard conditions
        global state1 = State("Air",state1.p,state9.mdot*state9.T+state10.mdot*state10.T,state1.mdot;phase = "gas",y_N2 = state9.mdot*state9.y_N2+state10.mdot*state10.y_N2,x_N2 = x_N2,liquid_fraction = state10.liquid_fraction)
        
        #check convergence
        println("\r",i);flush(stdout)
        #println(state1) 

        push!(solutions,state1)
        if length(solutions) > 1 && compare(solutions[end-1],solutions[end])
            break
        end
        if i == 20
            @error("Did not converge after 20 iterations")
        end
    end
    state1H = State("Essotherm650",stateOil_out1.p,(stateOil_out1.T*oil_distribution[1]+stateOil_out2.T*oil_distribution[2]),(stateOil_out1.mdot+stateOil_out2.mdot))
    return state1,state2A,state2B,state2C,state2,state4,state5,state6,state7,state9,state10,state1H,yield,T4,T9
end

function discharge_cycle(state1R::AirState,stateOil_in::OilState,oil_distribution,propane_max,methanol_max,pinch_coldbox,pinch_superheaters,pressure_after_pump,expander_pressures,η_e,η_pump,pressure_loss)
    number_of_expanders = 3
    if length(expander_pressures) != number_of_expanders && length(oil_distribution) != number_of_expanders
        error("The length of the 'expander_pressures' and the 'oil_distribution' must be equal to the number of expanders.")
    end

    state2R = isentropic_cryopump(state1R,pressure_after_pump,η_pump)
    state3R = heater_coldstorage(state2R,propane_max,pinch_coldbox,pressure_loss)
    state4R = heater_coldstorage(state3R,methanol_max,pinch_coldbox,pressure_loss)
    
    #distribution calculated from reference paper
    stateOil_in1 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[1])
    stateOil_in2 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[2])
    stateOil_in3 = State(stateOil_in.fluid,stateOil_in.p,stateOil_in.T,stateOil_in.mdot*oil_distribution[3])
        
    #calculate discharge process with guessed state 5R
    T_guess_5R = 430 #[K]
    p_5R = state4R.p-state4R.p*pressure_loss
    state_5R_guess = State("Air",p_5R,T_guess_5R,state4R.mdot;phase = "gas",y_N2 = state4R.y_N2,x_N2 = state4R.x_N2,liquid_fraction = state4R.liquid_fraction)

    #discharging
    state6R,stateOil_out1 = intercooler("Heat",state_5R_guess,stateOil_in1,pinch_superheaters,pressure_loss)
    state7R = isentropic_expander(state6R,expander_pressures[1],η_e)
    state8R,stateOil_out2 = intercooler("Heat",state7R,stateOil_in2,pinch_superheaters,pressure_loss)
    state9R = isentropic_expander(state8R,expander_pressures[2],η_e)
    state10R,stateOil_out3 = intercooler("Heat",state9R,stateOil_in3,pinch_superheaters,pressure_loss)
    state11R = isentropic_expander(state10R,expander_pressures[3],η_e)
    state12R = State("Air",state11R.p-state11R.p*pressure_loss,state4R.T+pinch_coldbox,state11R.mdot;phase = state11R.phase,y_N2 = state11R.y_N2,x_N2 = state11R.x_N2,liquid_fraction = state11R.liquid_fraction)
    
    #energy balance over heat exchanger
    h_5R = state11R.h-state12R.h+state4R.h
    T_5R = find_T_5R(h_5R,p_5R,state_5R_guess.y_N2)
        
    #define state5R
    state5R = State("Air",p_5R,ustrip(T_5R),state4R.mdot;phase = "gas",y_N2 = state4R.y_N2,x_N2 = state4R.x_N2,liquid_fraction = state4R.liquid_fraction)

    #discharging
    global state6R,stateOil_out1 = intercooler("Heat",state_5R_guess,stateOil_in1,pinch_superheaters,pressure_loss)
    global state7R = isentropic_expander(state6R,expander_pressures[1],η_e)
    global state8R,stateOil_out2 = intercooler("Heat",state7R,stateOil_in2,pinch_superheaters,pressure_loss)
    global state9R = isentropic_expander(state8R,expander_pressures[2],η_e)
    global state10R,stateOil_out3 = intercooler("Heat",state9R,stateOil_in3,pinch_superheaters,pressure_loss)
    global state11R = isentropic_expander(state10R,expander_pressures[3],η_e)
    global state12R = State("Air",state11R.p-state11R.p*pressure_loss,state4R.T+pinch_coldbox,state11R.mdot;phase = state11R.phase,y_N2 = state11R.y_N2,x_N2 = state11R.x_N2,liquid_fraction = state11R.liquid_fraction)

    state4H = State("Essotherm650",stateOil_out1.p,(stateOil_out1.T*oil_distribution[1]+stateOil_out2.T*oil_distribution[2]+stateOil_out3.T*oil_distribution[3]),(stateOil_out1.mdot+stateOil_out2.mdot+stateOil_out3.mdot))

    return state1R,state2R,state3R,state4R,state5R,state6R,state7R,state8R,state9R,state10R,state11R,state12R,state4H
end
