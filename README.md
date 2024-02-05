# StochasticStir.jl

Repository for the implementation of the stochastic stirring forcing and jet drag following
[Barnes and Hartmann, 2011](https://journals.ametsoc.org/view/journals/atsc/68/12/jas-d-11-039.1.xml?tab_body=fulltext-display#bib14)
and [Vallis et al. 2004](https://journals.ametsoc.org/view/journals/atsc/61/3/1520-0469_2004_061_0264_amasdm_2.0.co_2.xml)
for [SpeedyWeather.jl](https://github.com/SpeedyWeather/SpeedyWeather.jl).

## Background

See 
- https://github.com/SpeedyWeather/SpeedyWeather.jl/issues/394
- Issue https://github.com/SpeedyWeather/StochasticStir.jl/issues/3

## Installation

Via Julia's package manager (opened with `]`)

```julia
julia>]add https://github.com/SpeedyWeather/StochasticStir.jl
```
or via `using Pkg; Pkg.add("https://github.com/SpeedyWeather/StochasticStir.jl")`.
This package is not registered in Julia's general registry, so you have to specify
the whole URL as above. Requires SpeedyWeather.jl v0.8 or higher.

## Usage

This package exports `StochasticStirring` and `JetDrag`. The respective terms can be added
to SpeedyWeather.jl's `BarotropicModel` or `ShallowWaterModel` as follows. Feel free to
raise an issue here for any related questions.

```julia
using SpeedyWeather
using StochasticStir

# model components
spectral_grid = SpectralGrid(trunc=42,nlev=1)   
forcing = StochasticStirring(spectral_grid,latitude=45,strength=7e-11)
drag = JetDrag(spectral_grid,time_scale=Day(6))
initial_conditions = StartFromRest()

# construct the model and initialize
model = BarotropicModel(;spectral_grid,initial_conditions,forcing,drag)
simulation = initialize!(model)

# now run and store output
run!(simulation,period=Day(100),output=true)
```

and you can get a list of options for `JetDrag` and `StochasticStirring` by typing
```julia
julia>?JetDrag
julia>?StochasticStirring
```
where `?` opens the help.

## Gallery

A visualisation of the above simulation 

https://github.com/SpeedyWeather/StochasticStir.jl/assets/25530332/3924a9fd-8515-4628-8128-f555e6054616
