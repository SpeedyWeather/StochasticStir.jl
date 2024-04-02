"""
Defines the JetDrag drag term for SpeedyWeather that applies
a zonal drag to the BarotropicModel or ShallowWaterModel following
a jet at specified latitude, strength and width with time scale time_scale.
$(TYPEDFIELDS)"""
Base.@kwdef struct JetDrag{NF} <: SpeedyWeather.AbstractDrag

    # DIMENSIONS from SpectralGrid
    "Spectral resolution as max degree of spherical harmonics"
    trunc::Int

    # OPTIONS
    "Relaxation time scale τ"
    time_scale::Second = Day(6)

    "Jet strength [m/s]"
    u₀::NF = 20

    "latitude of Gaussian jet [˚N]"
    latitude::NF = 30

    "Width of Gaussian jet [˚]"
    width::NF = 6

    # TO BE INITIALISED
    "Relaxation back to reference vorticity"
    ζ₀::LTM{Complex{NF}} = zeros(LTM{Complex{NF}}, trunc+2, trunc+1)
end

"""
$(TYPEDSIGNATURES)
Generator function of a JetDrag struct taking resolution and number format from 
a SpectralGrid. Other options are provided as keyword arguments."""
function JetDrag(SG::SpectralGrid; kwargs...)
    return JetDrag{SG.NF}(; SG.trunc, kwargs...)
end

"""
$(TYPEDSIGNATURES)
Extends the initialize! function for a JetDrag. Precomputes the zonal velocity
of the jet, transforms to spectral and takes the curl to get the corresponding
relative vorticity of the jet, which is used to apply the drag to the vorticity
equation."""
function SpeedyWeather.initialize!( drag::JetDrag,
                                    model::ModelSetup)

    (;spectral_grid, geometry) = model
    (;Grid, NF, nlat_half) = spectral_grid
    lat = geometry.latds        # latitude in ˚N for every grid point

    # zonal velocity of the jet, allocate following grid and number format NF of model
    u = zeros(Grid{NF}, nlat_half)

    # Gaussian in latitude, calculated for every grid point ij
    for ij in eachindex(u)
        u[ij] = drag.u₀ * exp(-(lat[ij]-drag.latitude)^2/(2*drag.width^2))
    end

    # to spectral space, reusing the precomputed spectral transform from model
    û = SpeedyTransforms.spectral(u, model.spectral_transform)
    v̂ = zero(û)
    SpeedyTransforms.curl!(drag.ζ₀, û, v̂, model.spectral_transform)
    return nothing
end

# function barrier to unpack only what's needed
function SpeedyWeather.drag!( 
    diagn::DiagnosticVariablesLayer,
    progn::PrognosticVariablesLayer,
    drag::JetDrag,
    time::DateTime,
    model::ModelSetup,
)
    SpeedyWeather.drag!(diagn, progn, drag, model.geometry)
end

"""
$(TYPEDSIGNATURES)
Extends the drag! function for JetDrag. Adds the -r*(ζ-ζ₀) term
to the vorticity tendency.
"""
function SpeedyWeather.drag!(
    diagn::DiagnosticVariablesLayer,
    progn::PrognosticVariablesLayer,
    drag::JetDrag,
    geometry::Geometry,
)
    (; vor) = progn
    (; vor_tend) = diagn.tendencies
    (; ζ₀) = drag
    
    # 1/time_scale [1/s] but scaled with radius as is the vorticity equation
    (; radius) = geometry
    r = radius/drag.time_scale.value

    # add the -r*(ζ-ζ₀) term to the vorticity tendency that already includes the forcing term
    @inbounds for lm in eachindex(vor, vor_tend, ζ₀)
        vor_tend[lm] -= r*(vor[lm] - ζ₀[lm])
    end
end