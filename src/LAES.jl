module LAES
# Packages
using CoolProp
using Unitful: m, J, kPa, Pa, kg, K, W, °C, l, bar, MW, g, uconvert, °, ustrip, s, kW, hr,mol,kJ
using Plots
using DataFrames
using CurveFit

# load in diffent function files
include("extra_file.jl")
include("components.jl")

export isentropic_compressor
export intercooler
export isentropic_cryoexpander
export separator
export coldbox
export isentropic_cryopump
export heater_coldstorage
export isentropic_expander
end
