# StochasticStir.jl

Repository for the implementation of the stochastic stirring forcing and jet drag following
[Barnes and Hartmann, 2011](https://journals.ametsoc.org/view/journals/atsc/68/12/jas-d-11-039.1.xml?tab_body=fulltext-display#bib14)
and [Vallis et al. 2004](https://journals.ametsoc.org/view/journals/atsc/61/3/1520-0469_2004_061_0264_amasdm_2.0.co_2.xml).

See also https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/394

To be used like

```julia
using SpeedyWeather
using StochasticStirr

# model components
spectral_grid = SpectralGrid(trunc=42,nlev=1)   
forcing = StochasticStirring(spectral_grid,latitude=45)
drag = JetDrag(spectral,time_scale=Day(6))
initial_conditions = StartFromRest()

# construct the model and initialize
model = BarotropicModel(;spectral_grid,initial_conditions,forcing,drag)
simulation = initialize!(model)

# now run
run!(simulation)
```
