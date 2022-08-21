using LAES
using Test
using Unitful: ustrip
#########################################################################
# Use the reference paper of Guizzi to validate the component functions
#########################################################################

# Data charging cycle
const state1_guizzi = State("Air",100000,296.24,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
const state2A_guizzi = State("Air",1480000,687.74,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
const state2B_guizzi = State("Air",1465000,308.15,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
const state2C_guizzi = State("Air",18098000,682.00,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
const state2_guizzi = State("Air",17917000,308.15,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
const state3_guizzi = State("Air",17738000,245.80,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
const state4_guizzi = State("Air",17561000,98.00,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
const state5_guizzi = State("Air",102000,78.91,1;phase = "2phase",y_N2 = 0.93,x_N2 = 0.77,liquid_fraction = 0.842)
const state6_guizzi = State("Air",102000,78.91,0.842;phase = "liquid",y_N2 = 0,x_N2 = 0.77,liquid_fraction = 0.842)
const state7_guizzi = State("Air",102000,78.91,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
const state8_guizzi = State("Air",101000,237.80,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
const state9_guizzi = State("Air",100000,286.28,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
const state10_guizzi = State("Air",100000,298.15,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)

# Data discharge cycle
const state1R_guizzi = State("Air",100000,78.74,0.842;phase = "liquid",y_N2 = 0,x_N2 = 0.77,liquid_fraction = 1)
const state2R_guizzi = State("Air",6500000,81.89,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state3R_guizzi = State("Air",6435000,209.00,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state4R_guizzi = State("Air",6371000,283.00,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state5R_guizzi = State("Air",6307000,436.27,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state6R_guizzi = State("Air",6244000,616.42,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state7R_guizzi = State("Air",1590000,450.55,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state8R_guizzi = State("Air",1574000,616.42,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state9R_guizzi = State("Air",401000,451.23,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state10R_guizzi = State("Air",397000,616.42,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state11R_guizzi = State("Air",101000,451.42,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
const state12R_guizzi = State("Air",100000,288.00,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)

# Data thermal oil (Essotherm650)
const state1H_guizzi = State("Essotherm650",100000,626.42,0.999)
const state2H_guizzi = State("Essotherm650",100000,288.15,0.999)
const state3H_guizzi = State("Essotherm650",100000,626.42,0.999)
const state4H_guizzi = State("Essotherm650",100000,460.71,0.999)

    # Thermal oil distribution charging (calculated with the enthalpies of the reference paper + applying energy conservation)
oil1 = (707.454-308.727)/(0.999*(849.94-26.95))
oil2 = (705.204-281.71)/(0.999*(849.94-26.95))
oil_distribution = [oil1 oil2]

    # Thermal oil distribution discharge
oil1R = (628.96-436.35)/(0.999*(849.94-395.31))
oil2R = (628.96-454.68)/(0.999*(849.94-395.31))
oil3R = (629.01-456.27)/(0.999*(849.94-395.31))
        # sum isn't equal to 1 (probably an error)!!!!! --> rescale
oil1R_new = oil1R/sum([oil1R,oil2R,oil3R])  
oil2R_new = oil2R/sum([oil1R,oil2R,oil3R])  
oil3R_new = oil3R/sum([oil1R,oil2R,oil3R])  
oil_distribution_R = [oil1R_new oil2R_new oil3R_new]

# Data coolants 
const state1CS_guizzi = State("Propane",100000,93,1.019)
const state2CS_guizzi = State("Propane",100000,214,1.019)
const state3CS_guizzi = State("Methanol",100000,214,0.437)
const state4CS_guizzi = State("Methanol",100000,288,0.437)

#########################################################################
# Tests charging cycle
#########################################################################

# Compressors
state2A = isentropic_compressor(state1_guizzi,1480000,0.85)
T_difference_2A = abs((state2A.T -state2A_guizzi.T)/state2A_guizzi.T)

state2C = isentropic_compressor(state2B_guizzi,18098000,0.85)
T_difference_2C = abs((state2C.T -state2C_guizzi.T)/state2C_guizzi.T)

@testset "Compressor" begin
    @test T_difference_2A < 0.01
    @test T_difference_2C < 0.01
end

# Intercoolers
stateOil_in1 = State(state2H_guizzi.fluid,state2H_guizzi.p,state2H_guizzi.T,state2H_guizzi.mdot*oil_distribution[1])
stateOil_in2 = State(state2H_guizzi.fluid,state2H_guizzi.p,state2H_guizzi.T,state2H_guizzi.mdot*oil_distribution[2])
    
state2B,stateOil_out1 = intercooler("Cool",state2A_guizzi,stateOil_in1,10,0.01)
state2,stateOil_out2 = intercooler("Cool",state2C_guizzi,stateOil_in2,10,0.01)

    # Outlet temperature of air isn't 308.15K !!!! --> pinch becomes 20K
state2B,stateOil_out1 = intercooler("Cool",state2A_guizzi,stateOil_in1,20,0.01)
T_difference_2B = abs((state2B.T -state2B_guizzi.T)/state2B_guizzi.T)
p_difference_2B = abs((state2B.p -state2B_guizzi.p)/state2B_guizzi.p)

state2,stateOil_out2 = intercooler("Cool",state2C_guizzi,stateOil_in2,20,0.01)  
T_difference_2 = abs((state2.T -state2_guizzi.T)/state2_guizzi.T)
p_difference_2 = abs((state2.p -state2_guizzi.p)/state2_guizzi.p)

state1H = State("Essotherm650",stateOil_out1.p,(stateOil_out1.T*oil_distribution[1]+stateOil_out2.T*oil_distribution[2]),(stateOil_out1.mdot+stateOil_out2.mdot))
T_difference_1H = abs((state1H.T -state1H_guizzi.T)/state1H_guizzi.T)


stateOil_in1_R = State(state1H_guizzi.fluid,state1H_guizzi.p,state1H_guizzi.T,state1H_guizzi.mdot*oil_distribution_R[1])
stateOil_in2_R = State(state1H_guizzi.fluid,state1H_guizzi.p,state1H_guizzi.T,state1H_guizzi.mdot*oil_distribution_R[2])
stateOil_in3_R = State(state1H_guizzi.fluid,state1H_guizzi.p,state1H_guizzi.T,state1H_guizzi.mdot*oil_distribution_R[3])
        
state6R,stateOil_out1_R = intercooler("Heat",state5R_guizzi,stateOil_in1_R,10,0.01)
T_difference_6R = abs((state6R.T -state6R_guizzi.T)/state6R_guizzi.T)
p_difference_6R = abs((state6R.p -state6R_guizzi.p)/state6R_guizzi.p)
state8R,stateOil_out2_R = intercooler("Heat",state7R_guizzi,stateOil_in2_R,10,0.01)
T_difference_8R = abs((state8R.T -state8R_guizzi.T)/state8R_guizzi.T)
p_difference_8R = abs((state8R.p -state8R_guizzi.p)/state8R_guizzi.p)
state10R,stateOil_out3_R = intercooler("Heat",state9R_guizzi,stateOil_in3_R,10,0.01)
T_difference_10R = abs((state10R.T -state10R_guizzi.T)/state10R_guizzi.T)
p_difference_10R = abs((state10R.p -state10R_guizzi.p)/state10R_guizzi.p)

state4H = State("Essotherm650",stateOil_out1_R.p,(stateOil_out1_R.T*oil_distribution_R[1]+stateOil_out2_R.T*oil_distribution_R[2]+stateOil_out3_R.T*oil_distribution_R[3]),(stateOil_out1_R.mdot+stateOil_out2_R.mdot+stateOil_out3_R.mdot))
T_difference_4H = abs((state4H.T -state4H_guizzi.T)/state4H_guizzi.T)
mdot_difference_4H = abs((state4H.mdot -state4H_guizzi.mdot)/state4H_guizzi.mdot)

@testset "Intercoolers" begin
    @test T_difference_2B < 0.01
    @test p_difference_2B < 0.01 #rounding error
    @test T_difference_2 < 0.01
    @test p_difference_2 < 0.01 #rounding error
    @test T_difference_1H < 0.01
    @test T_difference_6R < 0.01
    @test p_difference_6R < 0.01 #rounding error
    @test T_difference_8R < 0.01
    @test p_difference_8R < 0.01 #rounding error
    @test T_difference_10R < 0.01
    @test p_difference_10R < 0.01 #rounding error
    @test T_difference_4H < 0.01
    @test mdot_difference_4H < 0.01
end

# Cryo-expander
state5 =  isentropic_cryoexpander(state4_guizzi,102000,0.7)
T_difference_5 = abs((state5.T -state5_guizzi.T)/state5_guizzi.T)
y_N2_difference_5 = abs((state5.y_N2 -state5_guizzi.y_N2)/state5_guizzi.y_N2)
x_N2_difference_5 = abs((state5.x_N2 -state5_guizzi.x_N2)/state5_guizzi.x_N2)
liquid_fraction_difference_5 = abs((state5.liquid_fraction -state5_guizzi.liquid_fraction)/state5_guizzi.liquid_fraction)

@testset "Cryo-expander" begin
    @test T_difference_5 < 0.01
    @test y_N2_difference_5 < 0.01 
    @test x_N2_difference_5 < 0.01
    @test liquid_fraction_difference_5 < 0.01
end

# Separator
state6,state7 = separator(state5_guizzi)

@testset "Separator" begin
    @test state6.p == state6_guizzi.p
    @test state7.p == state7_guizzi.p
    @test state6.T == state6_guizzi.T
    @test state7.T == state7_guizzi.T
    @test state6.phase == state6_guizzi.phase
    @test state7.phase == state7_guizzi.phase
    @test state6.y_N2 == state6_guizzi.y_N2
    @test state7.y_N2 == state7_guizzi.y_N2
    @test state6.x_N2 == state6_guizzi.x_N2
    @test state7.x_N2 == state7_guizzi.x_N2
    @test state6.liquid_fraction == state6_guizzi.liquid_fraction
    @test state7.liquid_fraction == state7_guizzi.liquid_fraction
end

#########################################################################
# Tests storage
#########################################################################

state1R = storage_tank(state6_guizzi,100000)
T_difference_1R = abs((state1R.T -state1R_guizzi.T)/state1R_guizzi.T)

@testset "Storage" begin
    @test T_difference_1R < 0.01
    @test state1R.mdot == state1R_guizzi.mdot
    @test state1R.x_N2 == state1R_guizzi.x_N2
    @test state1R.y_N2 == state1R_guizzi.y_N2
    @test state1R.liquid_fraction == state1R_guizzi.liquid_fraction
end

#########################################################################
# Tests discharge cycle
#########################################################################

# Cryo-pump
state2R = isentropic_cryopump(state1R_guizzi,6500000,0.7)
T_difference_2R = abs((state2R.T -state2R_guizzi.T)/state2R_guizzi.T)

@testset "Cryo-pump" begin
    @test T_difference_2R < 0.01
    @test state2R.mdot == state2R_guizzi.mdot
    @test state2R.x_N2 == state2R_guizzi.x_N2
    @test state2R.y_N2 == state2R_guizzi.y_N2
    @test state2R.liquid_fraction == state2R_guizzi.liquid_fraction
end

# Heaters
state3R = heater_coldstorage(state2R_guizzi,state2CS_guizzi,5,0.01)
state4R = heater_coldstorage(state3R_guizzi,state4CS_guizzi,5,0.01)

@testset "Heaters" begin
    @test state3R.T == state3R_guizzi.T
    @test state4R.T == state4R_guizzi.T
    @test state3R.phase == state3R_guizzi.phase
    @test state4R.phase == state4R_guizzi.phase
    @test state3R.y_N2 == state3R_guizzi.y_N2
    @test state4R.y_N2 == state4R_guizzi.y_N2
    @test state3R.x_N2 == state3R_guizzi.x_N2
    @test state4R.x_N2 == state4R_guizzi.x_N2
    @test state3R.liquid_fraction == state3R_guizzi.liquid_fraction
    @test state4R.liquid_fraction == state4R_guizzi.liquid_fraction
end

# Expanders
state7R = isentropic_expander(state6R_guizzi,state7R_guizzi.p,0.85)
T_difference_7R = abs((state7R.T -state7R_guizzi.T)/state7R_guizzi.T)

state9R = isentropic_expander(state8R_guizzi,state9R_guizzi.p,0.85)
T_difference_9R = abs((state9R.T -state9R_guizzi.T)/state9R_guizzi.T)

state11R = isentropic_expander(state10R_guizzi,state11R_guizzi.p,0.85)
T_difference_11R = abs((state11R.T -state11R_guizzi.T)/state11R_guizzi.T)


@testset "Expanders" begin
    @test T_difference_7R < 0.01
    @test T_difference_9R < 0.01
    @test T_difference_11R < 0.01
    @test state7R.phase == state7R_guizzi.phase
    @test state9R.phase == state9R_guizzi.phase
    @test state11R.phase == state11R_guizzi.phase
    @test state7R.y_N2 == state7R_guizzi.y_N2
    @test state9R.y_N2 == state9R_guizzi.y_N2
    @test state11R.y_N2 == state11R_guizzi.y_N2
    @test state7R.x_N2 == state7R_guizzi.x_N2
    @test state9R.x_N2 == state9R_guizzi.x_N2
    @test state11R.x_N2 == state11R_guizzi.x_N2
    @test state7R.liquid_fraction == state7R_guizzi.liquid_fraction
    @test state9R.liquid_fraction == state9R_guizzi.liquid_fraction
    @test state11R.liquid_fraction == state11R_guizzi.liquid_fraction
end

#########################################################################
# Tests extra functions
#########################################################################

T_5R = find_T_5R(state5R_guizzi.h,state5R_guizzi.p,state5R_guizzi.y_N2)
T_difference_5R = abs((ustrip(T_5R) -state5R_guizzi.T)/state5R_guizzi.T)

@testset "Extras" begin
    @test T_difference_5R < 0.01
end