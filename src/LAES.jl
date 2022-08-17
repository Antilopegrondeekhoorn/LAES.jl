module LAES
# Packages
using CoolProp
using Unitful: m, J, kPa, Pa, kg, K, W, °C, l, MW, g, uconvert, °, ustrip, s, kW, hr,mol,kJ
using DataFrames
using CurveFit

# load in diffent function files
include("components.jl")
include("helper_functions.jl")
include("simulation.jl")

export State
export isentropic_T
export T_real
export h_Essotherm650
export T_Essotherm650
export interpolate_K
export flash_calculation
export T_compressed_air
export pinch_coldbox_optimal
export pinch_coldbox_p_less_optimal
export pinch_coldbox_p_more_optimal
export compare
export find_T_5R

export storage_tank

export isentropic_compressor
export intercooler
export isentropic_cryoexpander
export separator
export coldbox
export isentropic_cryopump
export heater_coldstorage
export isentropic_expander

export charging_cycle
export discharge_cycle
end
