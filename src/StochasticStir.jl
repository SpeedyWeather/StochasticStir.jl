module StochasticStir

using SpeedyWeather

export StochasticStirring, JetDrag

include("stochastic_stirring.jl")
include("jet_drag.jl")

end
