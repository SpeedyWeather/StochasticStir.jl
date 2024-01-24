module StochasticStir

using SpeedyWeather
using DocStringExtensions

export StochasticStirring, JetDrag

include("stochastic_stirring.jl")
include("jet_drag.jl")

end
