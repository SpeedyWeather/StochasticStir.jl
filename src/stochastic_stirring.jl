export StochasticStirring
Base.@kwdef struct StochasticStirring{NF} <: SpeedyWeather.AbstractForcing{NF}
    
    # DIMENSIONS from SpectralGrid
    "Spectral resolution as max degree of spherical harmonics"
    trunc::Int
    
    "Number of latitude rings, used for latitudinal mask"
    nlat::Int

    
    # OPTIONS
    "Decorrelation time scale τ [days]"
    decorrelation_time::Second = Day(2)

    "Stirring strength A [1/s²]"
    strength::Float64 = 1e-11

    "Stirring latitude [˚N]"
    latitude::Float64 = 45

    "Stirring width [˚]"
    width::Float64 = 24

    "Minimum degree of spherical harmonics to force"
    lmin::Int = 8

    "Maximum degree of spherical harmonics to force"
    lmax::Int = 40

    "Minimum order of spherical harmonics to force"
    mmin::Int = 4

    "Maximum order of spherical harmonics to force"
    mmax::Int = lmax
    
    # TO BE INITIALISED
    "Stochastic stirring term S"
    S::LowerTriangularMatrix{Complex{NF}} = zeros(LowerTriangularMatrix{Complex{NF}},trunc+2,trunc+1)
    
    "a = A*sqrt(1 - exp(-2dt/τ)), the noise factor times the stirring strength [1/s²]"
    a::Base.RefValue{NF} = Ref(zero(NF))
        
    "b = exp(-dt/τ), the auto-regressive factor [1]"
    b::Base.RefValue{NF} = Ref(zero(NF))
        
    "Latitudinal mask, confined to mid-latitude storm track by default [1]"
    lat_mask::Vector{NF} = zeros(NF,nlat)
end

function StochasticStirring(SG::SpectralGrid;kwargs...)
    (;trunc,Grid,nlat_half) = SG
    nlat = RingGrids.get_nlat(Grid,nlat_half)
    return StochasticStirring{SG.NF}(;trunc,nlat,kwargs...)
end

function SpeedyWeather.initialize!( forcing::StochasticStirring,
                                    model::ModelSetup)
    
    # precompute forcing strength, scale with radius^2 as is the vorticity equation
    (;radius) = model.spectral_grid
    A = radius^2 * forcing.strength
    
    # precompute noise and auto-regressive factor, packed in RefValue for mutability
    dt = model.time_stepping.Δt_sec
    τ = forcing.decorrelation_time.value    # in seconds
    forcing.a[] = A*sqrt(1 - exp(-2dt/τ))
    forcing.b[] = exp(-dt/τ)
    
    # precompute the latitudinal mask
    (;Grid,nlat_half) = model.spectral_grid
    latd = RingGrids.get_latd(Grid,nlat_half)
    
    for j in eachindex(forcing.lat_mask)
        # Gaussian centred at forcing.latitude of width forcing.width
        forcing.lat_mask[j] = exp(-(forcing.latitude-latd[j])^2/forcing.width^2*2)
    end

    return nothing
end

function SpeedyWeather.forcing!(diagn::DiagnosticVariablesLayer,
                                progn::PrognosticVariablesLayer,
                                forcing::StochasticStirring,
                                time::DateTime,
                                model::ModelSetup)
    SpeedyWeather.forcing!(diagn,forcing,model.spectral_transform)
end

function SpeedyWeather.forcing!(diagn::DiagnosticVariablesLayer,
                                forcing::StochasticStirring{NF},
                                spectral_transform::SpectralTransform) where NF
    
    # noise and auto-regressive factors
    a = forcing.a[]    # = sqrt(1 - exp(-2dt/τ))
    b = forcing.b[]    # = exp(-dt/τ)
    
    (;S) = forcing
    lmax,mmax = size(S)
    @inbounds for m in 1:mmax
        for l in m:lmax
            if (forcing.mmin <= m <= forcing.mmax) &&
                (forcing.lmin <= l <= forcing.lmax)
                # Barnes and Hartmann, 2011 Eq. 2
                Qi = 2rand(Complex{NF}) - (1 + im)   # ~ [-1,1] in complex
                S[l,m] = a*Qi + b*S[l,m]
            end
        end
    end

    # to grid-point space
    S_grid = diagn.dynamics_variables.a_grid
    SpeedyTransforms.gridded!(S_grid,S,spectral_transform)
    
    # mask everything but mid-latitudes
    RingGrids._scale_lat!(S_grid,forcing.lat_mask)
    
    # back to spectral space
    (;vor_tend) = diagn.tendencies
    SpeedyTransforms.spectral!(vor_tend,S_grid,spectral_transform)
    SpeedyTransforms.spectral_truncation!(vor_tend)    # set lmax+1 to zero
    
    return nothing
end