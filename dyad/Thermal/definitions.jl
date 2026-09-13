# Julia helper functions visible to Dyad components in this library.
# Loaded automatically from `dyad/definitions.jl` by the generated module.jl.

"""
    sun_height_deg(t, latitude_rad, longitude_rad, day_of_year, time_zone)

Solar altitude angle above horizon, in degrees, at simulation time `t` (seconds since
midnight UTC of `day_of_year`).

Combines the equation-of-time-free model from ShipSIM's `SolarIrradiation`:
- Declination δ = 23.45° · sin(2π · (284 + N) / 365), where N = day_of_year.
- True mean time `tmt = (t/3600) - time_zone + 12·longitude_rad/π`.
- Hour angle ω = (tmt - 12) · 15° = π·(tmt - 12)/12 rad.
- sin(altitude) = sin(δ)·sin(φ) + cos(δ)·cos(φ)·cos(ω)
"""
function sun_height_deg(t, latitude_rad, longitude_rad, day_of_year, time_zone)
    delta = deg2rad(23.45) * sin(2 * pi * (284 + day_of_year) / 365)
    tmt = t / 3600 - time_zone + 12 * longitude_rad / pi
    omega = pi * (tmt - 12) / 12
    sin_alt = sin(delta) * sin(latitude_rad) +
              cos(delta) * cos(latitude_rad) * cos(omega)
    return rad2deg(asin(sin_alt))
end

"""
    sun_vector_world(t, latitude_rad, longitude_rad, day_of_year, time_zone, irradiance_ref)

Sun direction vector in world axes (X = East, Y = North, Z = Up), magnitude scaled by
`irradiance_ref` so the vector encodes both direction and clear-sky irradiance.

Returns `[Sx, Sy, Sz]` pointing from sun → ground (so -S is the direction *to* the sun).
At zero altitude (horizon) the vector is horizontal; at zenith it's straight down (-Z).
"""
function sun_vector_world(t, latitude_rad, longitude_rad, day_of_year, time_zone, irradiance_ref)
    delta = deg2rad(23.45) * sin(2 * pi * (284 + day_of_year) / 365)
    tmt = t / 3600 - time_zone + 12 * longitude_rad / pi
    omega = pi * (tmt - 12) / 12
    sin_alt = sin(delta) * sin(latitude_rad) +
              cos(delta) * cos(latitude_rad) * cos(omega)
    cos_alt = sqrt(1 - sin_alt^2 + 1e-12)
    sin_az = -cos(delta) * sin(omega) / cos_alt
    cos_az = (sin(delta) * cos(latitude_rad) -
              cos(delta) * sin(latitude_rad) * cos(omega)) / cos_alt
    # Sun direction (from sun to ground), in (East=X, North=Y, Up=Z) world frame.
    # `sin_az`/`cos_az` describe the direction *to* the sun, so all three components are
    # negated: at solar noon in the northern hemisphere the sun stands due south and the
    # rays travel north and downward, giving Y > 0 and Z < 0.
    # Below-horizon zeroing is left to the Dyad side (`ifelse(SunHeight > 0, ...)`).
    return [-irradiance_ref * sin_az * cos_alt,
            -irradiance_ref * cos_az * cos_alt,
            -irradiance_ref * sin_alt]
end

"""
    discretize_cylinder_areas(R, N) -> Vector{Float64}

Equal-area radial discretization of a cylinder of radius `R` into `N` rings.
Returns the area of each ring (all equal: `π·R²/N`) — kept as a function so the user
can override later if a non-uniform discretization is desired.
"""
function discretize_cylinder_areas(R::Real, N::Integer)
    An = pi * R^2 / N
    return fill(An, N)
end

"""
    discretize_cylinder_conductances(R, L, k, N) -> Vector{Float64}

Equal-area radial discretization. Returns inter-ring conductances such that
`G[i] = Ht[i]·L·k` where `Ht[i] = 2π/log(Rm[i+1]/Rm[i])` for `i < N` and
`Ht[N] = 2π/log(Re[N]/Rm[N])` is the boundary conductance from the last ring to the
outer surface. Length of returned array is `N`.
"""
function discretize_cylinder_conductances(R::Real, L::Real, k::Real, N::Integer)
    An = pi * R^2 / N
    Re = zeros(N)
    Rm = zeros(N)
    Re[1] = sqrt(An / pi)
    Rm[1] = Re[1] / sqrt(2)
    for i in 2:N
        Re[i] = sqrt((An + pi * Re[i-1]^2) / pi)
        Rm[i] = (Re[i-1] + Re[i]) / 2
    end
    Ht = zeros(N)
    if N == 1
        Ht[1] = 1e30   # effectively-infinite (the upstream uses `Modelica.Constants.inf`).
    else
        for i in 1:(N-1)
            Ht[i] = 2 * pi / log(Rm[i+1] / Rm[i])
        end
        Ht[N] = 2 * pi / log(Re[N] / Rm[N])
    end
    return Ht .* (L * k)
end
