using SpeedyWeather
using StochasticStir
using Test

@testset "StochasticStir.jl" begin
    spectral_grid = SpectralGrid(trunc=31, nlev=1)
    
    drag = JetDrag(spectral_grid, time_scale=Day(6))
    forcing = StochasticStirring(spectral_grid)
    initial_conditions = StartFromRest()

    # with barotropic model
    model = BarotropicModel(;spectral_grid, initial_conditions, forcing, drag)
    simulation = initialize!(model)

    run!(simulation, period=Day(5))
    @test simulation.model.feedback.nars_detected == false

    # with shallow water model
    model = ShallowWaterModel(;spectral_grid, initial_conditions, forcing, drag)
    simulation = initialize!(model)

    run!(simulation, period=Day(5))
    @test simulation.model.feedback.nars_detected == false
end

@testset "Initialize at various resolutions" begin
    for trunc in (31, 42, 63, 85, 127, 170, 255, 340)
        spectral_grid = SpectralGrid(trunc=trunc, nlev=1)
        drag = JetDrag(spectral_grid, time_scale=Day(6))
        forcing = StochasticStirring(spectral_grid)
        initial_conditions = StartFromRest()
    
        # with barotropic model
        model = BarotropicModel(;spectral_grid, initial_conditions, forcing, drag)
        simulation = initialize!(model)
    end
end

@testset "Mask off" begin
    spectral_grid = SpectralGrid(trunc=31, nlev=1)
    drag = JetDrag(spectral_grid, time_scale=Day(6))
    
    forcing = StochasticStirring(spectral_grid, mask=false)
    model = BarotropicModel(;spectral_grid, initial_conditions, forcing, drag)
    simulation = initialize!(model)

    run!(simulation, period=Day(5))
    @test simulation.model.feedback.nars_detected == false
end