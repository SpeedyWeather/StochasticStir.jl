Base.@kwdef struct JetDrag{NF} <: SpeedyWeather.AbstractDrag{NF}

    # DIMENSIONS from SpectralGrid
    "Spectral resolution as max degree of spherical harmonics"
    trunc::Int

    # OPTIONS
    "Relaxation time scale τ"
    time_scale::Second = Day(6)

    "Jet strength [m/s]"
    u₀::Float64 = 20

    "latitude of Gaussian jet [˚N]"
    latitude::Float64 = 30

    "Width of Gaussian jet [˚]"
    width::Float64 = 6

    # TO BE INITIALISED
    "Relaxation back to reference vorticity"
    ζ₀::LowerTriangularMatrix{Complex{NF}} = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
end

function JetDrag(SG::SpectralGrid;kwargs...)
    return JetDrag{SG.NF}(;SG.trunc,kwargs...)
end

function SpeedyWeather.initialize!( drag::JetDrag,
                                    model::ModelSetup)

    (;spectral_grid, geometry) = model
    (;Grid,NF,nlat_half) = spectral_grid
    u = zeros(Grid{NF},nlat_half)

    lat = geometry.latds

    for ij in eachindex(u)
        u[ij] = drag.u₀ * exp(-(lat[ij]-drag.latitude)^2/(2*drag.width^2))
    end

    û = SpeedyTransforms.spectral(u,one_more_degree=true)
    v̂ = zero(û)
    SpeedyTransforms.curl!(drag.ζ₀,û,v̂,model.spectral_transform)
    return nothing
end

function SpeedyWeather.drag!(   diagn::DiagnosticVariablesLayer,
                                progn::PrognosticVariablesLayer,
                                drag::JetDrag,
                                time::DateTime,
                                model::ModelSetup)

    (;vor) = progn
    (;vor_tend) = diagn.tendencies
    (;ζ₀) = drag
    
    (;radius) = model.spectral_grid
    r = radius/drag.time_scale.value
    for lm in eachindex(vor,vor_tend,ζ₀)
        vor_tend[lm] -= r*(vor[lm] - ζ₀[lm])
    end
end