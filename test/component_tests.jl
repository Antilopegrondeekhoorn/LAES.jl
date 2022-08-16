#using helper_functions
#using components
using LAES
using Test

# Use the reference paper of Guizzi to validate the functions
state1_guizzi = State("Air",100000,296.24,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state2A_guizzi = State("Air",1480000,687.74,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state2B_guizzi = State("Air",1465000,308.15,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state2C_guizzi = State("Air",18098000,682.00,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state2_guizzi = State("Air",17917000,308.15,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state3_guizzi = State("Air",17738000,245.80,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state4_guizzi = State("Air",17561000,98.00,1;phase = "gas",y_N2 = 0.795,x_N2 = 0,liquid_fraction = 0)
state5_guizzi = State("Air",102000,78.91,1;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
state6_guizzi = State("Air",102000,78.91,0.842;phase = "liquid",y_N2 = 0,x_N2 = 0.77,liquid_fraction = 0.842)
state7_guizzi = State("Air",102000,78.91,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
state8_guizzi = State("Air",101000,237.80,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
state9_guizzi = State("Air",100000,286.28,0.158;phase = "gas",y_N2 = 0.93,x_N2 = 0,liquid_fraction = 0)
state10_guizzi = State("Air",100000,298.15,0.842;phase = "gas",y_N2 = 0.77,x_N2 = 0,liquid_fraction = 0)

# Compressors
state2A = isentropic_compressor(state1_guizzi,1480000,0.85)
T_difference_2A = abs((state2A.T -state2A_guizzi.T)/state2A_guizzi.T)

state2C = isentropic_compressor(state2B_guizzi,18098000,0.85)
T_difference_2C = abs((state2C.T -state2C_guizzi.T)/state2C_guizzi.T)

@testset "Compressor" begin
    @test T_difference_2A < 0.01
    @test T_difference_2C < 0.01
end
