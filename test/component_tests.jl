#using helper_functions
#using components
using LAES
using Test

#########################################################################
# Use the reference paper of Guizzi to validate the functions
#########################################################################

# Data charging cycle
state1_guizzi = State("Air",100000,296.24,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state2A_guizzi = State("Air",1480000,687.74,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state2B_guizzi = State("Air",1465000,308.15,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state2C_guizzi = State("Air",18098000,682.00,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state2_guizzi = State("Air",17917000,308.15,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state3_guizzi = State("Air",17738000,245.80,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state4_guizzi = State("Air",17561000,98.00,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state5_guizzi = State("Air",102000,78.91,1;phase = "2phase",y_N2 = 0.93,x_N2 = 0.77,liquid_fraction = 0.842)
state6_guizzi = State("Air",102000,78.91,0.842;phase = "liquid",y_N2 = 0,x_N2 = 0.77,liquid_fraction = 0.842)
state7_guizzi = State("Air",102000,78.91,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
state8_guizzi = State("Air",101000,237.80,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
state9_guizzi = State("Air",100000,286.28,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
state10_guizzi = State("Air",100000,298.15,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)

# Data discharge cycle
state1R_guizzi = State("Air",100000,78.74,0.842;phase = "liquid",y_N2 = 0,x_N2 = 0.77,liquid_fraction = 1)
state2R_guizzi = State("Air",6500000,81.89,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state3R_guizzi = State("Air",6435000,209.00,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state4R_guizzi = State("Air",6371000,283.00,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state5R_guizzi = State("Air",6307000,436.27,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state6R_guizzi = State("Air",6244000,616.42,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state7R_guizzi = State("Air",1590000,450.55,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state8R_guizzi = State("Air",1574000,616.42,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state9R_guizzi = State("Air",401000,451.23,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state10R_guizzi = State("Air",397000,616.42,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state11R_guizzi = State("Air",101000,451.42,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)
state12R_guizzi = State("Air",100000,288.00,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)

# Data thermal oil (Essotherm650)
state1H_guizzi = State("Essotherm650",100000,626.42,0.999)
state2H_guizzi = State("Essotherm650",100000,288.15,0.999)
state3H_guizzi = State("Essotherm650",100000,626.42,0.999)
state4H_guizzi = State("Essotherm650",100000,460.71,0.999)

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
state1C_guizzi = State("Propane",100000,93,1.019)
state2C_guizzi = State("Propane",100000,214,1.019)
state3C_guizzi = State("Methanol",100000,214,0.437)
state4C_guizzi = State("Methanol",100000,288,0.437)

#########################################################################
# Tests charging cycle
#########################################################################

# State functions


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

    # Outlet temperature of air isn't 308.15K (probably an error)!!!! --> pinch becomes 20K
state2B,stateOil_out1 = intercooler("Cool",state2A_guizzi,stateOil_in1,20,0.01)
T_difference_2B = abs((state2B.T -state2B_guizzi.T)/state2B_guizzi.T)
p_difference_2B = abs((state2B.p -state2B_guizzi.p)/state2B_guizzi.p)

state2,stateOil_out2 = intercooler("Cool",state2C_guizzi,stateOil_in2,20,0.01)  
T_difference_2 = abs((state2.T -state2_guizzi.T)/state2_guizzi.T)
p_difference_2 = abs((state2.p -state2_guizzi.p)/state2_guizzi.p)

state1H = State("Essotherm650",stateOil_out1.p,(stateOil_out1.T*oil_distribution[1]+stateOil_out2.T*oil_distribution[2]),(stateOil_out1.mdot+stateOil_out2.mdot))
T_difference_1H = abs((state1H.T -state1H_guizzi.T)/state1H_guizzi.T)

@testset "Intercoolers" begin
    @test T_difference_2B < 0.01
    @test p_difference_2B < 0.01 #rounding error
    @test T_difference_2 < 0.01
    @test p_difference_2 < 0.01 #rounding error
    @test T_difference_1H < 0.01
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

# Heaters
state3R = heater_coldstorage(state2R_guizzi,state2C_guizzi,5,0.01)
state4R = heater_coldstorage(state3R_guizzi,state4C_guizzi,5,0.01)
