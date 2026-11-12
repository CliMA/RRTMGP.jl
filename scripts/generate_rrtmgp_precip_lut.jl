"""
    generate_rrtmgp_precip_lut.jl

Generate standalone snow and rain optical property lookup tables
(`rrtmgp-precip-lw-bnd.nc` and `rrtmgp-precip-sw-bnd.nc`) for RRTMGP.jl based on:
  - Feng, J., Menzel, R., & Paynter, D. (2025).
    "Toward Transparency and Consistency: An Open-Source Optics Parameterization
    for Clouds and Precipitation." Journal of Advances in Modeling Earth Systems,
    17(3), e2024MS004478. https://doi.org/10.1029/2024MS004478
  - Yang, P., Bi, L., Baum, B. A., Liou, K.-N., Kattawar, G. W., Mishchenko, M. I.,
    & Cole, B. (2013). "Spectrally Consistent Scattering, Absorption, and
    Polarization Properties of Atmospheric Ice Crystals at Wavelengths from 0.2
    to 100 μm." Journal of the Atmospheric Sciences, 70(1), 330–347.

Usage:
  julia --threads=auto --project=scripts scripts/generate_rrtmgp_precip_lut.jl \\
      <input_data_dir> <output_dir> [snow_habit] [snow_roughness]

where `<input_data_dir>` contains:
  - `MIR/<snow_habit>/Rough<NNN>/isca.dat`
  - `FIR/<snow_habit>/Rough<NNN>/isca.dat`
  - `refractive_index` (Segelstein 1981 liquid water refractive index)
  - `solar_flux.csv` (solar irradiance spectrum)
"""

using NCDatasets
using Printf

const CM_TO_UM = 1.0e4
const M_TO_MICRON = 1.0e6
const RHO_ICE_G_M3 = 917.0e3   # g / m³ (Feng et al. 2025)
const RHO_LIQ_G_M3 = 996.0e3   # g / m³ (Feng et al. 2025)

# RRTMGP 16 Longwave bands [cm⁻¹]
const RRTMGP_LW_BANDS = [
    10.0 250.0
    250.0 500.0
    500.0 630.0
    630.0 700.0
    700.0 820.0
    820.0 980.0
    980.0 1080.0
    1080.0 1180.0
    1180.0 1390.0
    1390.0 1480.0
    1480.0 1800.0
    1800.0 2080.0
    2080.0 2250.0
    2250.0 2390.0
    2390.0 2680.0
    2680.0 3250.0
]

# RRTMGP 14 Shortwave bands [cm⁻¹]
const RRTMGP_SW_BANDS = [
    820.0 2680.0
    2680.0 3250.0
    3250.0 4000.0
    4000.0 4650.0
    4650.0 5150.0
    5150.0 6150.0
    6150.0 7700.0
    7700.0 8050.0
    8050.0 12850.0
    12850.0 16000.0
    16000.0 22650.0
    22650.0 29000.0
    29000.0 38000.0
    38000.0 50000.0
]

# ---------------------------------------------------------------------------
# 1. 1D Linear Interpolation with Boundary Clamping
# ---------------------------------------------------------------------------

"""
    interp1d_clamped_scalar(x::AbstractVector{Float64}, y::AbstractVector{Float64}, xq::Float64)

1D linear interpolation with constant clamping outside `[first(x), last(x)]`,
assuming `x` is strictly increasing.
"""
@inline function interp1d_clamped_scalar(
    x::AbstractVector{Float64},
    y::AbstractVector{Float64},
    xq::Float64,
)
    n = length(x)
    @assert n >= 2
    @inbounds if xq <= x[1]
        return y[1]
    elseif xq >= x[n]
        return y[n]
    else
        idx = clamp(searchsortedlast(x, xq), 1, n - 1)
        x0, x1 = x[idx], x[idx + 1]
        t = (xq - x0) / (x1 - x0)
        return muladd(t, y[idx + 1] - y[idx], y[idx])
    end
end

function interp1d_clamped(
    x::AbstractVector{Float64},
    y::AbstractVector{Float64},
    xq::AbstractVector{Float64},
)
    return [interp1d_clamped_scalar(x, y, xqi) for xqi in xq]
end

# ---------------------------------------------------------------------------
# 2. Spectral Utilities: Planck Function, Solar Spectrum, Refractive Index
# ---------------------------------------------------------------------------

"""
    planck_wn(wn::AbstractVector{Float64}, T::Float64)

Evaluate blackbody spectral intensity [W m⁻² sr⁻¹ cm] at wavenumber `wn` [cm⁻¹]
and temperature `T` [K].
"""
function planck_wn(wn::AbstractVector{Float64}, T::Float64)
    c1 = 1.191066e-8 # W / m² / sr / cm⁻⁴
    c2 = 1.438833    # K / cm⁻¹
    return @. c1 * wn^3 / (exp(c2 * wn / T) - 1.0)
end

"""
    read_solar_spectrum(path::AbstractString)

Read solar irradiance spectrum `(wavenum [cm⁻¹], solar_flux)` from CSV.
"""
function read_solar_spectrum(path::AbstractString)
    wavenum = Float64[]
    solar = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            isempty(strip(line)) && continue
            parts = split(line, ',')
            push!(wavenum, parse(Float64, parts[1]))
            push!(solar, parse(Float64, parts[2]))
        end
    end
    return wavenum, solar
end

"""
    read_refractive_index(wavelength_um::AbstractVector{Float64}, path::AbstractString)

Read Segelstein (1981) liquid water complex refractive index `m = n - i*k` and
interpolate onto `wavelength_um` [μm].
"""
function read_refractive_index(
    wavelength_um::AbstractVector{Float64},
    path::AbstractString,
)
    freq = Float64[]
    n_re = Float64[]
    k_im = Float64[]
    open(path, "r") do io
        for _ in 1:3
            readline(io)
        end
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            parts = split(s)
            push!(freq, parse(Float64, parts[1]))
            push!(n_re, parse(Float64, parts[2]))
            push!(k_im, parse(Float64, parts[3]))
        end
    end
    @assert issorted(freq)
    @assert minimum(wavelength_um) >= first(freq) &&
            maximum(wavelength_um) <= last(freq)
    n_interp = interp1d_clamped(freq, n_re, wavelength_um)
    k_interp = interp1d_clamped(freq, k_im, wavelength_um)
    return ComplexF64.(n_interp, .-k_interp)
end

"""
    create_size_list(var_min::Float64, var_max::Float64, n::Real = 10)

Quasi-logarithmic characteristic-size grid generator from Feng et al. (2025).
"""
function create_size_list(var_min::Float64, var_max::Float64, n::Real = 10)
    result = Float64[var_min]
    dr = var_min
    while result[end] < var_max
        if 2 * n * dr < result[end]
            dr = result[end] / n
        end
        push!(result, result[end] + dr)
    end
    return result
end

# ---------------------------------------------------------------------------
# 3. Wiscombe (1980) / Bohren & Huffman Mie Scattering for Homogeneous Spheres
# ---------------------------------------------------------------------------

function _lentz_dn(z::ComplexF64, N::Int)
    zinv = 2.0 / z
    alpha = (N + 0.5) * zinv
    aj = -(N + 1.5) * zinv
    alpha_j1 = aj + 1.0 / alpha
    alpha_j2 = aj
    ratio = alpha_j1 / alpha_j2
    runratio = alpha * ratio
    while abs(abs(ratio) - 1.0) > 1.0e-12
        aj = zinv - aj
        alpha_j1 = 1.0 / alpha_j1 + aj
        alpha_j2 = 1.0 / alpha_j2 + aj
        ratio = alpha_j1 / alpha_j2
        zinv = -zinv
        runratio = ratio * runratio
    end
    return -N / z + runratio
end

function _d_downwards!(D::Vector{ComplexF64}, z::ComplexF64, N::Int)
    last_D = _lentz_dn(z, N)
    @inbounds for n in N:-1:1
        last_D = n / z - 1.0 / (last_D + n / z)
        D[n] = last_D
    end
    return D
end

function _d_upwards!(D::Vector{ComplexF64}, z::ComplexF64, N::Int)
    exp_z = exp(-2im * z)
    D[1] = 0.0 + 0.0im
    if N >= 2
        D[2] =
            -1.0 / z + (1.0 - exp_z) / ((1.0 - exp_z) / z - 1im * (1.0 + exp_z))
        @inbounds for n in 3:N
            m_n = n - 1
            D[n] = 1.0 / (m_n / z - D[n - 1]) - m_n / z
        end
    end
    return D
end

function _small_mie(m::ComplexF64, x::Float64)
    m2 = m * m
    x2 = x * x
    D =
        m2 + 2.0 + (1.0 - 0.7 * m2) * x2 -
        (8.0 * m2 * m2 - 385.0 * m2 + 350.0) * x2 * x2 / 1400.0 +
        2im * (m2 - 1.0) * x^3 * (1.0 - 0.1 * x2) / 3.0
    ahat1 =
        2im * (m2 - 1.0) / 3.0 *
        (1.0 - 0.1 * x2 + (4.0 * m2 + 5.0) * x2 * x2 / 1400.0) / D
    bhat1 =
        1im * x2 * (m2 - 1.0) / 45.0 * (1.0 + (2.0 * m2 - 5.0) / 70.0 * x2) /
        (1.0 - (2.0 * m2 - 5.0) / 30.0 * x2)
    ahat2 =
        1im * x2 * (m2 - 1.0) / 15.0 * (1.0 - x2 / 14.0) /
        (2.0 * m2 + 3.0 - (2.0 * m2 - 7.0) / 14.0 * x2)
    T = abs2(ahat1) + abs2(bhat1) + (5.0 / 3.0) * abs2(ahat2)
    g = real(ahat1 * conj(ahat2 + bhat1)) / T
    qsca = 6.0 * x2 * x2 * T
    qext =
        imag(m) == 0.0 ? qsca :
        6.0 * x * real(ahat1 + bhat1 + 5.0 * ahat2 / 3.0)
    return qext, qsca, g
end

"""
    mie_scalar(m::ComplexF64, x::Float64, D_buf::Vector{ComplexF64})

Compute `(qext, qsca, g)` for a sphere with complex refractive index `m` and
size parameter `x`, following Wiscombe (1980) / `miepython`.
"""
function mie_scalar(m::ComplexF64, x::Float64, D_buf::Vector{ComplexF64})
    if real(m) > 0.0 && abs(m) * x < 0.1
        return _small_mie(m, x)
    end
    nstop = Int(floor(x + 4.05 * cbrt(x) + 2.0)) + 1
    N_D = nstop + 1
    length(D_buf) < N_D && resize!(D_buf, N_D)

    n_re = real(m)
    kappa = abs(imag(m))
    z = m * x
    if n_re < 1.0 ||
       n_re > 10.0 ||
       kappa > 10.0 ||
       x * kappa >= 3.9 - 10.8 * n_re + 13.78 * n_re^2
        _d_downwards!(D_buf, z, N_D)
    else
        _d_upwards!(D_buf, z, N_D)
    end

    psi_nm1 = sin(x)
    psi_n = psi_nm1 / x - cos(x)
    xi_nm1 = ComplexF64(psi_nm1, cos(x))
    xi_n = ComplexF64(psi_n, cos(x) / x + sin(x))

    sum_ext = 0.0
    sum_sca = 0.0
    sum_g = 0.0
    a_prev = 0.0 + 0.0im
    b_prev = 0.0 + 0.0im

    nmax = nstop - 1
    @inbounds for n in 1:nmax
        Dn = D_buf[n + 1]
        temp_a = Dn / m + n / x
        a_n = (temp_a * psi_n - psi_nm1) / (temp_a * xi_n - xi_nm1)
        temp_b = Dn * m + n / x
        b_n = (temp_b * psi_n - psi_nm1) / (temp_b * xi_n - xi_nm1)

        cn = 2.0 * n + 1.0
        sum_ext += cn * (real(a_n) + real(b_n))
        sum_sca += cn * (abs2(a_n) + abs2(b_n))

        c2n = cn / (n * (n + 1.0))
        asy2 = c2n * real(a_n * conj(b_n))
        if n > 1
            np = n - 1.0
            c1n = np * (np + 2.0) / (np + 1.0)
            asy1 = c1n * real(a_prev * conj(a_n) + b_prev * conj(b_n))
            sum_g += asy1 + asy2
        else
            sum_g += asy2
        end

        a_prev = a_n
        b_prev = b_n

        xi = (2 * n + 1) * xi_n / x - xi_nm1
        xi_nm1 = xi_n
        xi_n = xi
        psi_nm1 = psi_n
        psi_n = real(xi_n)
    end

    x2 = x * x
    qext = 2.0 * sum_ext / x2
    qsca = imag(m) == 0.0 ? qext : 2.0 * sum_sca / x2
    g = 4.0 * sum_g / (qsca * x2)
    return qext, qsca, g
end

# ---------------------------------------------------------------------------
# 4. Yang et al. (2013) Ice Single-Scattering Reader
# ---------------------------------------------------------------------------

function read_isca(path::AbstractString)
    wavelength = Float64[]
    d = Float64[]
    v = Float64[]
    s = Float64[]
    qext = Float64[]
    ssa = Float64[]
    asy = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            parts = split(line)
            length(parts) < 7 && continue
            # Columns in Yang et al. (2013) isca.dat:
            # 1: wavelength [μm], 2: D_max [μm], 3: volume [μm³],
            # 4: projected area [μm²], 5: Q_ext [-], 6: SSA [-], 7: g [-]
            push!(wavelength, parse(Float64, parts[1]))
            push!(d, parse(Float64, parts[2]))
            push!(v, parse(Float64, parts[3]))
            push!(s, parse(Float64, parts[4]))
            push!(qext, parse(Float64, parts[5]))
            push!(ssa, parse(Float64, parts[6]))
            push!(asy, parse(Float64, parts[7]))
        end
    end
    return (; wavelength, d, v, s, qext, ssa, asy)
end

"""
    read_yang_ice_library(data_dir::AbstractString; habit="solid_column", roughness=50)

Read, validate, and merge the MIR and FIR `isca.dat` tables from Yang et al. (2013).
Returns a NamedTuple `(; wavenum, r_eff, d, s, v, qext, ssa, ext_cs, sca_cs, asy)`
where 2D arrays have shape `(nwav, nsize)` sorted by increasing wavenumber [cm⁻¹]
and increasing maximum dimension `d` [μm].
"""
function read_yang_ice_library(
    data_dir::AbstractString;
    habit::AbstractString = "solid_column",
    roughness::Int = 50,
)
    rgh_str = @sprintf("Rough%03d", roughness)
    mir_path = joinpath(data_dir, "MIR", habit, rgh_str, "isca.dat")
    fir_path = joinpath(data_dir, "FIR", habit, rgh_str, "isca.dat")

    mir = read_isca(mir_path)
    fir = read_isca(fir_path)

    d_mir_uniq = sort(unique(mir.d))
    d_fir_uniq = sort(unique(fir.d))
    d_mir_uniq == d_fir_uniq || error(
        "MIR and FIR maximum-dimension grids do not match in $mir_path vs $fir_path",
    )

    d_table = d_mir_uniq
    wl_table = sort(unique(vcat(mir.wavelength, fir.wavelength)))
    nwl = length(wl_table)
    nsize = length(d_table)

    s_table = fill(NaN, nsize)
    v_table = fill(NaN, nsize)
    qext_wl = zeros(Float64, nwl, nsize)
    ssa_wl = zeros(Float64, nwl, nsize)
    asy_wl = zeros(Float64, nwl, nsize)
    filled = falses(nwl, nsize)

    for tbl in (mir, fir)
        for k in eachindex(tbl.wavelength)
            wl = tbl.wavelength[k]
            dk = tbl.d[k]
            iwl = searchsortedfirst(wl_table, wl)
            isz = searchsortedfirst(d_table, dk)
            @assert iwl <= nwl && wl_table[iwl] == wl
            @assert isz <= nsize && d_table[isz] == dk

            if isnan(s_table[isz])
                s_table[isz] = tbl.s[k]
                v_table[isz] = tbl.v[k]
            else
                isapprox(s_table[isz], tbl.s[k]; rtol = 1e-5) || error(
                    "Inconsistent projected area at D=$dk μm: $(s_table[isz]) vs $(tbl.s[k])",
                )
                isapprox(v_table[isz], tbl.v[k]; rtol = 1e-5) || error(
                    "Inconsistent crystal volume at D=$dk μm: $(v_table[isz]) vs $(tbl.v[k])",
                )
            end

            # If a wavelength is duplicated across MIR and FIR, keep the first (MIR) entry
            if !filled[iwl, isz]
                qext_wl[iwl, isz] = tbl.qext[k]
                ssa_wl[iwl, isz] = tbl.ssa[k]
                asy_wl[iwl, isz] = tbl.asy[k]
                filled[iwl, isz] = true
            end
        end
    end

    all(filled) || error(
        "Incomplete (wavelength, size) grid when merging $mir_path and $fir_path: " *
        "$(count(!, filled)) missing entries out of $(length(filled))",
    )

    wavenum_unsorted = CM_TO_UM ./ wl_table
    iwav = sortperm(wavenum_unsorted)
    wavenum = wavenum_unsorted[iwav]

    qext = qext_wl[iwav, :]
    ssa = ssa_wl[iwav, :]
    asy = asy_wl[iwav, :]
    ext_cs = qext .* reshape(s_table, 1, nsize)
    sca_cs = ssa .* ext_cs
    r_eff = 0.75 .* v_table ./ s_table

    return (;
        wavenum,
        r_eff,
        d = d_table,
        s = s_table,
        v = v_table,
        qext,
        ssa,
        ext_cs,
        sca_cs,
        asy,
    )
end

# ---------------------------------------------------------------------------
# 5. Gamma PSD Integration and Spectral Band Averaging (Feng et al. 2025)
# ---------------------------------------------------------------------------

"""
    log_midpoint_quadrature(d_min::Float64, d_max::Float64, n_quad::Int)

Construct a log-spaced midpoint quadrature rule `(d_mid, dd)` on `[d_min, d_max]`
with `n_quad` sub-intervals, resolving both small particle sizes and wide
exponential tails accurately.
"""
function log_midpoint_quadrature(d_min::Float64, d_max::Float64, n_quad::Int)
    edges = exp.(range(log(d_min), log(d_max), length = n_quad + 1))
    d_mid = Vector{Float64}(undef, n_quad)
    dd = Vector{Float64}(undef, n_quad)
    @inbounds for j in 1:n_quad
        d_mid[j] = 0.5 * (edges[j] + edges[j + 1])
        dd[j] = edges[j + 1] - edges[j]
    end
    return d_mid, dd
end

"""
    gamma_psd_integrate(a, d_char, d_quad, dd_quad, v_quad, s_quad, ext_quad, sca_quad, asy_quad, rau)

Integrate single-particle optical properties over a Gamma PSD:
  `f(D) ∝ (D / b)^(a - 1) * exp(-D / b)` with `b = d_char / a`
using quadrature nodes `d_quad` and widths `dd_quad`.
Returns `(r_eff, s_mean, v_mean, ext, sca, ssa, asy)`.
"""
function gamma_psd_integrate(
    a::Float64,
    d_char::Float64,
    d_quad::Vector{Float64},
    dd_quad::Vector{Float64},
    v_quad::Vector{Float64},
    s_quad::Vector{Float64},
    ext_quad::Matrix{Float64}, # (nwav, nquad)
    sca_quad::Matrix{Float64}, # (nwav, nquad)
    asy_quad::Matrix{Float64}, # (nwav, nquad)
    rau::Float64,
)
    b = d_char / a
    nr = length(d_quad)
    nwav = size(ext_quad, 1)

    # Unnormalized Gamma weights w_j = (d_j / b)^(a - 1) * exp(-d_j / b) * dd_j.
    # The constant normalization factor (1 / (b * Gamma(a))) cancels out in all ratios.
    w = Vector{Float64}(undef, nr)
    int_f = 0.0
    int_fv = 0.0
    int_fs = 0.0
    @inbounds for j in 1:nr
        u = d_quad[j] / b
        wj = u^(a - 1.0) * exp(-u) * dd_quad[j]
        w[j] = wj
        int_f += wj
        int_fv += wj * v_quad[j]
        int_fs += wj * s_quad[j]
    end

    v_mean = int_fv / int_f
    s_mean = int_fs / int_f
    r_eff = 0.75 * v_mean / s_mean

    ext = zeros(Float64, nwav)
    sca = zeros(Float64, nwav)
    ssa = zeros(Float64, nwav)
    asy = zeros(Float64, nwav)

    norm_mass = M_TO_MICRON / (int_fv * rau) # m² / g
    @inbounds for iw in 1:nwav
        sum_ext = 0.0
        sum_sca = 0.0
        sum_asy_sca = 0.0
        for j in 1:nr
            wj = w[j]
            sc = max(sca_quad[iw, j], 0.0)
            as = max(asy_quad[iw, j], 0.0)
            sum_ext += wj * max(ext_quad[iw, j], 0.0)
            sum_sca += wj * sc
            sum_asy_sca += wj * as * sc
        end
        ext_iw = max(sum_ext * norm_mass, 0.0)
        sca_iw = clamp(sum_sca * norm_mass, 0.0, ext_iw)
        ssa_iw = ext_iw > 0.0 ? clamp(sca_iw / ext_iw, 0.0, 1.0) : 0.0
        asy_iw = sum_sca > 0.0 ? clamp(sum_asy_sca / sum_sca, 0.0, 1.0) : 0.0

        ext[iw] = ext_iw
        sca[iw] = sca_iw
        ssa[iw] = ssa_iw
        asy[iw] = asy_iw
    end

    return r_eff, s_mean, v_mean, ext, sca, ssa, asy
end

"""
    thin_average_bands(wavenum, source, band_limits, ext, ssa, asy)

Optically thin source-weighted band average (used for Longwave).
Weights extinction by `source`, scattering by `source * ext * ssa`, and
asymmetry parameter `g` by scattering optical depth × Planck source
(`source * ext * ssa`).
Input `ext`, `ssa`, `asy` have shape `(nr, nwav)`.
Returns `(ext_band, ssa_band, asy_band)` of shape `(nr, nband)`.
"""
function thin_average_bands(
    wavenum::AbstractVector{Float64},
    source::AbstractVector{Float64},
    band_limits::AbstractMatrix{Float64},
    ext::AbstractMatrix{Float64},
    ssa::AbstractMatrix{Float64},
    asy::AbstractMatrix{Float64},
)
    nr = size(ext, 1)
    nband = size(band_limits, 1)
    ext_b = zeros(Float64, nr, nband)
    ssa_b = zeros(Float64, nr, nband)
    asy_b = zeros(Float64, nr, nband)

    for ib in 1:nband
        wn_lwr, wn_upr = band_limits[ib, 1], band_limits[ib, 2]
        idx = findall(wn -> wn_lwr <= wn <= wn_upr, wavenum)
        @assert !isempty(idx)
        src_sum = sum(@view source[idx])
        @inbounds for ir in 1:nr
            e_sum = 0.0
            s_sum = 0.0
            g_sum = 0.0
            g_src_sum = 0.0
            for iw in idx
                sw = source[iw]
                e_iw = ext[ir, iw]
                sc_iw = e_iw * ssa[ir, iw]
                g_iw = asy[ir, iw]
                e_sum += sw * e_iw
                s_sum += sw * sc_iw
                g_sum += sw * sc_iw * g_iw
                g_src_sum += sw * g_iw
            end
            ext_val = max(e_sum / src_sum, 0.0)
            sca_val = clamp(s_sum / src_sum, 0.0, ext_val)
            ssa_val = ext_val > 0.0 ? clamp(sca_val / ext_val, 0.0, 1.0) : 0.0
            asy_val =
                s_sum > 0.0 ? clamp(g_sum / s_sum, 0.0, 1.0) :
                clamp(g_src_sum / src_sum, 0.0, 1.0)
            ext_b[ir, ib] = ext_val
            ssa_b[ir, ib] = ssa_val
            asy_b[ir, ib] = asy_val
        end
    end
    return ext_b, ssa_b, asy_b
end

"""
    thick_average_bands(wavenum, source, band_limits, ext, ssa, asy)

Edwards & Slingo (1996) optically thick source-weighted band average (used for Shortwave).
Input `ext`, `ssa`, `asy` have shape `(nr, nwav)`.
Returns `(ext_band, ssa_band, asy_band)` of shape `(nr, nband)`.
"""
function thick_average_bands(
    wavenum::AbstractVector{Float64},
    source::AbstractVector{Float64},
    band_limits::AbstractMatrix{Float64},
    ext::AbstractMatrix{Float64},
    ssa::AbstractMatrix{Float64},
    asy::AbstractMatrix{Float64},
)
    nr = size(ext, 1)
    nband = size(band_limits, 1)
    ext_b = zeros(Float64, nr, nband)
    ssa_b = zeros(Float64, nr, nband)
    asy_b = zeros(Float64, nr, nband)

    for ib in 1:nband
        wn_lwr, wn_upr = band_limits[ib, 1], band_limits[ib, 2]
        idx = findall(wn -> wn_lwr <= wn <= wn_upr, wavenum)
        @assert !isempty(idx)
        src_sum = sum(@view source[idx])
        @inbounds for ir in 1:nr
            e_sum = 0.0
            R_sum = 0.0
            g_sum = 0.0
            for iw in idx
                sw = source[iw]
                w0 = ssa[ir, iw]
                g0 = asy[ir, iw]
                s_val =
                    (w0 * g0 == 1.0) ? 0.0 :
                    sqrt(max(0.0, (1.0 - w0) / (1.0 - w0 * g0)))
                R_val = (1.0 - s_val) / (1.0 + s_val)
                e_sum += sw * ext[ir, iw]
                R_sum += sw * R_val
                g_sum += sw * g0
            end
            ext_val = max(e_sum / src_sum, 0.0)
            R_ave = R_sum / src_sum
            asy_val = clamp(g_sum / src_sum, 0.0, 1.0)
            ssa_val = if R_ave == 0.0
                0.0
            else
                clamp(
                    4.0 * R_ave / ((1.0 + R_ave)^2 - asy_val * (1.0 - R_ave)^2),
                    0.0,
                    1.0,
                )
            end
            ext_b[ir, ib] = ext_val
            ssa_b[ir, ib] = ssa_val
            asy_b[ir, ib] = asy_val
        end
    end
    return ext_b, ssa_b, asy_b
end

# ---------------------------------------------------------------------------
# 6. Full Snow (Yang 2013) and Rain (Mie) Band Table Generators
# ---------------------------------------------------------------------------

"""
    _fit_log_slope(x::AbstractVector{Float64}, y::AbstractVector{Float64})

Ordinary least-squares power-law exponent `β` fitting `log(y) ≈ α + β * log(x)`.
"""
function _fit_log_slope(x::AbstractVector{Float64}, y::AbstractVector{Float64})
    lx = log.(x)
    ly = log.(y)
    mx = sum(lx) / length(lx)
    my = sum(ly) / length(ly)
    num = sum((lx .- mx) .* (ly .- my))
    den = sum((lx .- mx) .^ 2)
    return num / den
end

"""
    compute_snow_band_tables(
        data_dir;
        habit = "solid_column",
        roughness = 50,
        a_psd = 1.0,
        lw_planck_temp = 233.0,
        n_quad = 2000,
    )

Compute LW and SW band tables `(; r_eff, ext_lw, ssa_lw, asy_lw, ext_sw, ssa_sw, asy_sw)`
for snow particles using Yang et al. (2013) single-scattering databases and
Feng et al. (2025) Gamma PSD + band averaging.
"""
function compute_snow_band_tables(
    data_dir::AbstractString;
    habit::AbstractString = "solid_column",
    roughness::Int = 50,
    a_psd::Float64 = 1.0,
    lw_planck_temp::Float64 = 233.0,
    n_quad::Int = 2000,
)
    yang = read_yang_ice_library(data_dir; habit, roughness)
    wavenum_in = yang.wavenum
    d_in = yang.d
    s_in = yang.s
    v_in = yang.v
    ext_cs_in = yang.ext_cs
    sca_cs_in = yang.sca_cs
    asy_in = yang.asy
    qext_in = yang.qext
    ssa_in = yang.ssa

    d_min_tbl = first(d_in)
    d_max_tbl = last(d_in)
    v_min_tbl, v_max_tbl = first(v_in), last(v_in)
    s_min_tbl, s_max_tbl = first(s_in), last(s_in)

    # Fit habit-specific large-crystal power-law exponents V ∝ D^β_V and S ∝ D^β_S
    # across the upper geometric-optics regime (D >= 0.7 * D_max) to continue V(D)
    # and S(D) smoothly beyond D_max = 10,000 μm without forcing spherical exponents.
    tail_idx = findall(>=(0.7 * d_max_tbl), d_in)
    beta_v = _fit_log_slope(@view(d_in[tail_idx]), @view(v_in[tail_idx]))
    beta_s = _fit_log_slope(@view(d_in[tail_idx]), @view(s_in[tail_idx]))

    cbrt_v_in = cbrt.(v_in)
    sqrt_s_in = sqrt.(s_in)
    asy_sca_in = asy_in .* sca_cs_in

    r_grid = create_size_list(0.5, 20000.0, 20)
    nr = length(r_grid)
    nwav_in = length(wavenum_in)

    r_out = zeros(Float64, nr)
    ext_inres = zeros(Float64, nr, nwav_in)
    ssa_inres = zeros(Float64, nr, nwav_in)
    asy_inres = zeros(Float64, nr, nwav_in)

    u_min = 1.0e-5 * max(1.0, a_psd)
    u_max = max(25.0, a_psd + 20.0 + 6.0 * sqrt(a_psd))

    Threads.@threads for i in 1:nr
        d_char = 2.0 * r_grid[i]
        b = d_char / a_psd
        d_quad, dd_quad = log_midpoint_quadrature(u_min * b, u_max * b, n_quad)

        v_quad = Vector{Float64}(undef, n_quad)
        s_quad = Vector{Float64}(undef, n_quad)
        @inbounds for j in 1:n_quad
            dj = d_quad[j]
            if dj < d_min_tbl
                ratio = dj / d_min_tbl
                v_quad[j] = v_min_tbl * ratio^3
                s_quad[j] = s_min_tbl * ratio^2
            elseif dj <= d_max_tbl
                v_quad[j] = interp1d_clamped_scalar(d_in, cbrt_v_in, dj)^3
                s_quad[j] = interp1d_clamped_scalar(d_in, sqrt_s_in, dj)^2
            else
                ratio = dj / d_max_tbl
                v_quad[j] = v_max_tbl * ratio^beta_v
                s_quad[j] = s_max_tbl * ratio^beta_s
            end
        end

        ext_quad = zeros(Float64, nwav_in, n_quad)
        sca_quad = zeros(Float64, nwav_in, n_quad)
        asy_quad = zeros(Float64, nwav_in, n_quad)
        for iw in 1:nwav_in
            ev = @view ext_cs_in[iw, :]
            sv = @view sca_cs_in[iw, :]
            asv = @view asy_sca_in[iw, :]
            qe_min, qe_max = qext_in[iw, 1], qext_in[iw, end]
            w0_min, w0_max = ssa_in[iw, 1], ssa_in[iw, end]
            g_min, g_max = asy_in[iw, 1], asy_in[iw, end]

            @inbounds for j in 1:n_quad
                dj = d_quad[j]
                sh = s_quad[j]
                if dj < d_min_tbl
                    e_val = qe_min * sh
                    s_val = w0_min * e_val
                    as_val = g_min
                elseif dj <= d_max_tbl
                    e_val = max(interp1d_clamped_scalar(s_in, ev, sh), 0.0)
                    s_val =
                        clamp(interp1d_clamped_scalar(s_in, sv, sh), 0.0, e_val)
                    as_val =
                        s_val > 0.0 ?
                        clamp(
                            interp1d_clamped_scalar(s_in, asv, sh) / s_val,
                            0.0,
                            1.0,
                        ) : 0.0
                else
                    # Geometric-optics asymptote for D > D_max: single-particle
                    # efficiencies (Q_ext, SSA, g) remain at their D_max limit.
                    e_val = qe_max * sh
                    s_val = w0_max * e_val
                    as_val = g_max
                end
                ext_quad[iw, j] = e_val
                sca_quad[iw, j] = s_val
                asy_quad[iw, j] = as_val
            end
        end

        r_eff, _, _, ext_i, _, ssa_i, asy_i = gamma_psd_integrate(
            a_psd,
            d_char,
            d_quad,
            dd_quad,
            v_quad,
            s_quad,
            ext_quad,
            sca_quad,
            asy_quad,
            RHO_ICE_G_M3,
        )
        r_out[i] = r_eff
        ext_inres[i, :] .= ext_i
        ssa_inres[i, :] .= ssa_i
        asy_inres[i, :] .= asy_i
    end

    @assert issorted(r_out)

    # Interpolate onto LW midpoint wavenumber grid (10.5:1.0:3249.5 cm⁻¹) and Planck-average
    wn_lw = collect(10.5:1.0:3249.5)
    src_lw = planck_wn(wn_lw, lw_planck_temp)
    ext_lw_hres, ssa_lw_hres, asy_lw_hres = _interp_optics_to_wavenum(
        wavenum_in,
        ext_inres,
        ssa_inres,
        asy_inres,
        wn_lw,
    )
    ext_lw, ssa_lw, asy_lw = thin_average_bands(
        wn_lw,
        src_lw,
        RRTMGP_LW_BANDS,
        ext_lw_hres,
        ssa_lw_hres,
        asy_lw_hres,
    )

    # Interpolate onto SW midpoint wavenumber grid (825.0:10.0:49995.0 cm⁻¹) and thick-average
    wn_sw = collect(825.0:10.0:49995.0)
    wn_sol, sol_flux = read_solar_spectrum(joinpath(data_dir, "solar_flux.csv"))
    src_sw = interp1d_clamped(wn_sol, sol_flux, wn_sw)
    ext_sw_hres, ssa_sw_hres, asy_sw_hres = _interp_optics_to_wavenum(
        wavenum_in,
        ext_inres,
        ssa_inres,
        asy_inres,
        wn_sw,
    )
    ext_sw, ssa_sw, asy_sw = thick_average_bands(
        wn_sw,
        src_sw,
        RRTMGP_SW_BANDS,
        ext_sw_hres,
        ssa_sw_hres,
        asy_sw_hres,
    )

    return (; r_eff = r_out, ext_lw, ssa_lw, asy_lw, ext_sw, ssa_sw, asy_sw)
end

function _interp_optics_to_wavenum(
    wn_in::Vector{Float64},
    ext_in::Matrix{Float64},
    ssa_in::Matrix{Float64},
    asy_in::Matrix{Float64},
    wn_out::Vector{Float64},
)
    nr = size(ext_in, 1)
    nw_out = length(wn_out)
    ext_out = zeros(Float64, nr, nw_out)
    ssa_out = zeros(Float64, nr, nw_out)
    asy_out = zeros(Float64, nr, nw_out)

    for i in 1:nr
        e_in = @view ext_in[i, :]
        es_in = e_in .* (@view ssa_in[i, :])
        esg_in = es_in .* (@view asy_in[i, :])
        @inbounds for j in 1:nw_out
            wn = wn_out[j]
            e_raw = interp1d_clamped_scalar(wn_in, e_in, wn)
            es_raw = interp1d_clamped_scalar(wn_in, es_in, wn)
            esg_raw = interp1d_clamped_scalar(wn_in, esg_in, wn)
            s_raw = e_raw > 0.0 ? es_raw / e_raw : 0.0
            g_raw = es_raw > 0.0 ? esg_raw / es_raw : 0.0
            ext_out[i, j] = max(e_raw, 0.0)
            ssa_out[i, j] = clamp(s_raw, 0.0, 1.0)
            asy_out[i, j] = clamp(g_raw, 0.0, 1.0)
        end
    end
    return ext_out, ssa_out, asy_out
end

"""
    compute_rain_band_tables(
        data_dir;
        a_psd = 1.0,
        lw_planck_temp = 250.0,
        dwn_lw = 1.0,
        dwn_sw = 10.0,
        n_quad = 2000,
    )

Compute LW and SW band tables `(; r_eff, ext_lw, ssa_lw, asy_lw, ext_sw, ssa_sw, asy_sw)`
for rain drops using exact Wiscombe (1980) Mie scattering with Segelstein (1981)
liquid water refractive indices and Feng et al. (2025) Gamma PSD + band averaging.
"""
function compute_rain_band_tables(
    data_dir::AbstractString;
    a_psd::Float64 = 1.0,
    lw_planck_temp::Float64 = 250.0,
    dwn_lw::Float64 = 1.0,
    dwn_sw::Float64 = 10.0,
    n_quad::Int = 2000,
)
    # For spheres under a Gamma PSD, r_eff = ((a + 2) / a) * r_char.
    # Choose characteristic radii so r_eff spans ~[2.4, 6000] μm, strictly bracketing [5, 5000] μm.
    r_char_min = 2.4 * a_psd / (a_psd + 2.0)
    r_char_max = 6000.0 * a_psd / (a_psd + 2.0)
    r_grid = create_size_list(r_char_min, r_char_max, 12)

    wn_lw = collect((10.0 + 0.5 * dwn_lw):dwn_lw:(3250.0 - 0.5 * dwn_lw))
    src_lw = planck_wn(wn_lw, lw_planck_temp)
    m_lw = read_refractive_index(
        CM_TO_UM ./ wn_lw,
        joinpath(data_dir, "refractive_index"),
    )

    wn_sw = collect((820.0 + 0.5 * dwn_sw):dwn_sw:(50000.0 - 0.5 * dwn_sw))
    wn_sol, sol_flux = read_solar_spectrum(joinpath(data_dir, "solar_flux.csv"))
    src_sw = interp1d_clamped(wn_sol, sol_flux, wn_sw)
    m_sw = read_refractive_index(
        CM_TO_UM ./ wn_sw,
        joinpath(data_dir, "refractive_index"),
    )

    r_out_lw, ext_lw_hres, ssa_lw_hres, asy_lw_hres =
        _compute_mie_psd_grid(r_grid, a_psd, wn_lw, m_lw; n_quad)
    r_out_sw, ext_sw_hres, ssa_sw_hres, asy_sw_hres =
        _compute_mie_psd_grid(r_grid, a_psd, wn_sw, m_sw; n_quad)
    @assert isapprox(r_out_lw, r_out_sw; rtol = 1e-12)

    ext_lw, ssa_lw, asy_lw = thin_average_bands(
        wn_lw,
        src_lw,
        RRTMGP_LW_BANDS,
        ext_lw_hres,
        ssa_lw_hres,
        asy_lw_hres,
    )
    ext_sw, ssa_sw, asy_sw = thick_average_bands(
        wn_sw,
        src_sw,
        RRTMGP_SW_BANDS,
        ext_sw_hres,
        ssa_sw_hres,
        asy_sw_hres,
    )

    return (; r_eff = r_out_lw, ext_lw, ssa_lw, asy_lw, ext_sw, ssa_sw, asy_sw)
end

function _compute_mie_psd_grid(
    r_grid::Vector{Float64},
    a_psd::Float64,
    wavenum::Vector{Float64},
    m_vec::Vector{ComplexF64};
    n_quad::Int = 2000,
)
    u_min = 1.0e-5 * max(1.0, a_psd)
    u_max = max(20.0, a_psd + 16.0 + 5.0 * sqrt(a_psd))

    r_min_single = min(0.01, 0.5 * u_min * (2.0 * minimum(r_grid) / a_psd))
    r_max_single = 0.5 * u_max * (2.0 * maximum(r_grid) / a_psd)
    r_split = min(2000.0, 0.5 * r_max_single)

    # Two-tier log-spaced single-particle droplet radius grid:
    # - 500 points across [r_min_single, r_split] to resolve Mie resonances densely
    # - 80 points across (r_split, r_max_single] in the smooth geometric-optics tail
    r_tier1 = exp.(range(log(r_min_single), log(r_split), length = 500))
    r_tier2 = exp.(range(log(r_split), log(r_max_single), length = 81))[2:end]
    r_in = vcat(r_tier1, r_tier2)
    n_single = length(r_in)
    nwav = length(wavenum)

    qext_in = zeros(Float64, nwav, n_single)
    qsca_in = zeros(Float64, nwav, n_single)
    asy_in = zeros(Float64, nwav, n_single)

    nmed = 1.0003
    Threads.@threads for iw in 1:nwav
        D_buf = Vector{ComplexF64}(undef, 512)
        m_w = m_vec[iw]
        wn = wavenum[iw]
        @inbounds for j in 1:n_single
            x = 2.0 * π * r_in[j] * wn / CM_TO_UM * nmed
            qe, qs, g = mie_scalar(m_w, x, D_buf)
            qext_in[iw, j] = max(qe, 0.0)
            qsca_in[iw, j] = clamp(qs, 0.0, max(qe, 0.0))
            asy_in[iw, j] = clamp(g, 0.0, 1.0)
        end
    end

    nr = length(r_grid)
    r_out = zeros(Float64, nr)
    ext_out = zeros(Float64, nr, nwav)
    ssa_out = zeros(Float64, nr, nwav)
    asy_out = zeros(Float64, nr, nwav)

    log_r_in = log.(r_in)
    qasy_in = asy_in .* qsca_in

    Threads.@threads for i in 1:nr
        d_char = 2.0 * r_grid[i]
        b = d_char / a_psd
        d_quad, dd_quad = log_midpoint_quadrature(u_min * b, u_max * b, n_quad)
        r_quad = 0.5 .* d_quad
        log_r_quad = log.(r_quad)

        v_quad = @. (4.0 / 3.0) * π * r_quad^3
        s_quad = @. π * r_quad^2

        ext_quad = zeros(Float64, nwav, n_quad)
        sca_quad = zeros(Float64, nwav, n_quad)
        asy_quad = zeros(Float64, nwav, n_quad)

        for iw in 1:nwav
            qe_row = @view qext_in[iw, :]
            qs_row = @view qsca_in[iw, :]
            qg_row = @view qasy_in[iw, :]
            @inbounds for j in 1:n_quad
                lr = log_r_quad[j]
                sh = s_quad[j]
                qe_val = max(interp1d_clamped_scalar(log_r_in, qe_row, lr), 0.0)
                qs_val = clamp(
                    interp1d_clamped_scalar(log_r_in, qs_row, lr),
                    0.0,
                    qe_val,
                )
                qg_val = max(interp1d_clamped_scalar(log_r_in, qg_row, lr), 0.0)
                ext_quad[iw, j] = qe_val * sh
                sca_quad[iw, j] = qs_val * sh
                asy_quad[iw, j] =
                    qs_val > 0.0 ? clamp(qg_val / qs_val, 0.0, 1.0) : 0.0
            end
        end

        r_eff, _, _, ext_i, _, ssa_i, asy_i = gamma_psd_integrate(
            a_psd,
            d_char,
            d_quad,
            dd_quad,
            v_quad,
            s_quad,
            ext_quad,
            sca_quad,
            asy_quad,
            RHO_LIQ_G_M3,
        )
        # Analytical self-check for spherical drops under a Gamma PSD:
        # r_eff_exact = ((a + 2) / (2 * a)) * d_char
        r_eff_exact = ((a_psd + 2.0) / (2.0 * a_psd)) * d_char
        isapprox(r_eff, r_eff_exact; rtol = 1e-4) || error(
            "Gamma PSD quadrature inaccurate at d_char=$d_char μm: got r_eff=$r_eff, expected $r_eff_exact",
        )

        r_out[i] = r_eff
        ext_out[i, :] .= ext_i
        ssa_out[i, :] .= ssa_i
        asy_out[i, :] .= asy_i
    end

    @assert issorted(r_out)
    return r_out, ext_out, ssa_out, asy_out
end

# ---------------------------------------------------------------------------
# 7. Resample onto Log-Uniform Effective Radius Grids for RRTMGP.jl
# ---------------------------------------------------------------------------

"""
    resample_to_log_grid(r_src, ext_src, ssa_src, asy_src, r_lwr, r_upr, nsize)

Resample band tables `(nr_src, nband)` onto a log-uniform effective radius grid:
  `r_grid[k] = r_lwr * (r_upr / r_lwr)^((k - 1) / (nsize - 1))` for `k = 1:nsize`.
Interpolates in `log(r_eff)` so that `compute_lookup_precip_props` in RRTMGP.jl
can perform O(1) direct index lookup in `log(r_eff)`.
Returns `(r_grid, ext_lut, ssa_lut, asy_lut)` with tables of shape `(nsize, nband)`.
"""
function resample_to_log_grid(
    r_src::Vector{Float64},
    ext_src::Matrix{Float64},
    ssa_src::Matrix{Float64},
    asy_src::Matrix{Float64},
    r_lwr::Float64,
    r_upr::Float64,
    nsize::Int,
)
    @assert issorted(r_src)
    first(r_src) <= r_lwr && last(r_src) >= r_upr || error(
        "Target LUT radius bounds [$r_lwr, $r_upr] μm exceed source PSD r_eff range " *
        "[$(first(r_src)), $(last(r_src))] μm",
    )

    log_r_src = log.(r_src)
    log_r_grid = range(log(r_lwr), log(r_upr), length = nsize)
    r_grid = exp.(collect(log_r_grid))
    r_grid[1] = r_lwr
    r_grid[end] = r_upr

    nband = size(ext_src, 2)
    ext_lut = zeros(Float32, nsize, nband)
    ssa_lut = zeros(Float32, nsize, nband)
    asy_lut = zeros(Float32, nsize, nband)

    for ib in 1:nband
        e_col = @view ext_src[:, ib]
        s_col = @view ssa_src[:, ib]
        g_col = @view asy_src[:, ib]
        for k in 1:nsize
            lr = log_r_grid[k]
            ext_lut[k, ib] =
                Float32(max(interp1d_clamped_scalar(log_r_src, e_col, lr), 0.0))
            ssa_lut[k, ib] = Float32(
                clamp(interp1d_clamped_scalar(log_r_src, s_col, lr), 0.0, 1.0),
            )
            asy_lut[k, ib] = Float32(
                clamp(interp1d_clamped_scalar(log_r_src, g_col, lr), 0.0, 1.0),
            )
        end
    end

    return Float32.(r_grid), ext_lut, ssa_lut, asy_lut
end

# ---------------------------------------------------------------------------
# 8. Write Standalone RRTMGP Precipitation NetCDF Files
# ---------------------------------------------------------------------------

"""
    write_precip_netcdf(
        dst_nc_path, λ_str, band_limits,
        r_sno, ext_sno, ssa_sno, asy_sno,
        r_rai, ext_rai, ssa_rai, asy_rai;
        snow_habit, snow_roughness, snow_a_psd, snow_lw_temp,
        rain_a_psd, rain_lw_temp,
    )

Write a standalone RRTMGP precipitation lookup NetCDF file (`rrtmgp-precip-lw-bnd.nc`
or `rrtmgp-precip-sw-bnd.nc`) containing `bnd_limits_wavenumber` plus log-uniform
snow and rain optical lookup tables.
"""
function write_precip_netcdf(
    dst_nc_path::AbstractString,
    λ_str::AbstractString,
    band_limits::AbstractMatrix{Float64}, # (nband, 2)
    r_sno::Vector{Float32},
    ext_sno::Matrix{Float32}, # (nsize_sno, nband)
    ssa_sno::Matrix{Float32},
    asy_sno::Matrix{Float32},
    r_rai::Vector{Float32},
    ext_rai::Matrix{Float32}, # (nsize_rai, nband)
    ssa_rai::Matrix{Float32},
    asy_rai::Matrix{Float32};
    snow_habit::AbstractString,
    snow_roughness::Int,
    snow_a_psd::Float64,
    snow_lw_temp::Float64,
    rain_a_psd::Float64,
    rain_lw_temp::Float64,
)
    mkpath(dirname(dst_nc_path))
    isfile(dst_nc_path) && rm(dst_nc_path)

    nband = size(band_limits, 1)
    bnd_lims_wn = Float32.(permutedims(band_limits)) # (2, nband) = ("pair", "nband")

    NCDataset(dst_nc_path, "c") do ds
        ds.attrib["title"] = "RRTMGP $(λ_str) precipitation (snow and rain) band optical lookup tables"
        ds.attrib["references"] = "Feng et al. (2025, JAMES, doi:10.1029/2024MS004478); Yang et al. (2013, JAS, doi:10.1175/JAS-D-12-039.1); Wiscombe (1980, Appl. Opt.); Segelstein (1981)"
        ds.attrib["generator"] = "RRTMGP.jl scripts/generate_rrtmgp_precip_lut.jl"
        ds.attrib["snow_habit"] = String(snow_habit)
        ds.attrib["snow_roughness"] = @sprintf("Rough%03d", snow_roughness)
        ds.attrib["snow_a_psd"] = snow_a_psd
        ds.attrib["rain_a_psd"] = rain_a_psd
        if λ_str == "longwave"
            ds.attrib["snow_lw_planck_temp_K"] = snow_lw_temp
            ds.attrib["rain_lw_planck_temp_K"] = rain_lw_temp
        end

        defDim(ds, "nband", nband)
        defDim(ds, "pair", 2)
        defDim(ds, "nsize_sno", length(r_sno))
        defDim(ds, "nsize_rai", length(r_rai))

        defVar(
            ds,
            "bnd_limits_wavenumber",
            bnd_lims_wn,
            ("pair", "nband"),
            attrib = Dict(
                "standard_name" => "bnd_limits_wavenumber",
                "description" => "Band wavenumber limits",
                "units" => "cm-1",
            ),
        )

        defVar(
            ds,
            "radsno_lwr",
            Float32(first(r_sno)),
            (),
            attrib = Dict(
                "standard_name" => "radsno_lwr",
                "description" => "Snow particle effective radius lower bound for band LUT log-interpolation",
                "units" => "microns",
            ),
        )
        defVar(
            ds,
            "radsno_upr",
            Float32(last(r_sno)),
            (),
            attrib = Dict(
                "standard_name" => "radsno_upr",
                "description" => "Snow particle effective radius upper bound for band LUT log-interpolation",
                "units" => "microns",
            ),
        )
        defVar(
            ds,
            "radsno",
            r_sno,
            ("nsize_sno",),
            attrib = Dict(
                "standard_name" => "radsno",
                "description" => "Snow particle effective radius log-uniform grid points",
                "units" => "microns",
            ),
        )
        defVar(
            ds,
            "extsno",
            ext_sno,
            ("nsize_sno", "nband"),
            attrib = Dict(
                "standard_name" => "extsno",
                "description" => "Band LUT $(λ_str) snow mass extinction coefficient (Feng et al. 2025 / Yang et al. 2013)",
                "units" => "m2/g",
            ),
        )
        defVar(
            ds,
            "ssasno",
            ssa_sno,
            ("nsize_sno", "nband"),
            attrib = Dict(
                "standard_name" => "ssasno",
                "description" => "Band LUT $(λ_str) snow single scattering albedo (Feng et al. 2025 / Yang et al. 2013)",
                "units" => "unitless",
            ),
        )
        defVar(
            ds,
            "asysno",
            asy_sno,
            ("nsize_sno", "nband"),
            attrib = Dict(
                "standard_name" => "asysno",
                "description" => "Band LUT $(λ_str) snow asymmetry parameter (Feng et al. 2025 / Yang et al. 2013)",
                "units" => "unitless",
            ),
        )

        defVar(
            ds,
            "radrai_lwr",
            Float32(first(r_rai)),
            (),
            attrib = Dict(
                "standard_name" => "radrai_lwr",
                "description" => "Rain particle effective radius lower bound for band LUT log-interpolation",
                "units" => "microns",
            ),
        )
        defVar(
            ds,
            "radrai_upr",
            Float32(last(r_rai)),
            (),
            attrib = Dict(
                "standard_name" => "radrai_upr",
                "description" => "Rain particle effective radius upper bound for band LUT log-interpolation",
                "units" => "microns",
            ),
        )
        defVar(
            ds,
            "radrai",
            r_rai,
            ("nsize_rai",),
            attrib = Dict(
                "standard_name" => "radrai",
                "description" => "Rain particle effective radius log-uniform grid points",
                "units" => "microns",
            ),
        )
        defVar(
            ds,
            "extrai",
            ext_rai,
            ("nsize_rai", "nband"),
            attrib = Dict(
                "standard_name" => "extrai",
                "description" => "Band LUT $(λ_str) rain mass extinction coefficient (Feng et al. 2025 Mie)",
                "units" => "m2/g",
            ),
        )
        defVar(
            ds,
            "ssarai",
            ssa_rai,
            ("nsize_rai", "nband"),
            attrib = Dict(
                "standard_name" => "ssarai",
                "description" => "Band LUT $(λ_str) rain single scattering albedo (Feng et al. 2025 Mie)",
                "units" => "unitless",
            ),
        )
        defVar(
            ds,
            "asyrai",
            asy_rai,
            ("nsize_rai", "nband"),
            attrib = Dict(
                "standard_name" => "asyrai",
                "description" => "Band LUT $(λ_str) rain asymmetry parameter (Feng et al. 2025 Mie)",
                "units" => "unitless",
            ),
        )
    end
    return dst_nc_path
end

"""
    generate_rrtmgp_precip_luts(;
        data_dir,
        lw_dst_nc,
        sw_dst_nc,
        snow_habit = "solid_column",
        snow_roughness = 50,
        snow_a_psd = 1.0,
        snow_lw_temp = 233.0,
        radsno_lwr = 5.0,
        radsno_upr = 1100.0,
        nsize_sno = 64,
        rain_a_psd = 1.0,
        rain_lw_temp = 250.0,
        radrai_lwr = 5.0,
        radrai_upr = 5000.0,
        nsize_rai = 64,
    )

End-to-end generator for standalone RRTMGP longwave and shortwave precipitation
lookup tables (`rrtmgp-precip-lw-bnd.nc` and `rrtmgp-precip-sw-bnd.nc`).
"""
function generate_rrtmgp_precip_luts(;
    data_dir::AbstractString,
    lw_dst_nc::AbstractString,
    sw_dst_nc::AbstractString,
    snow_habit::AbstractString = "solid_column",
    snow_roughness::Int = 50,
    snow_a_psd::Float64 = 1.0,
    snow_lw_temp::Float64 = 233.0,
    radsno_lwr::Float64 = 5.0,
    radsno_upr::Float64 = 1100.0,
    nsize_sno::Int = 64,
    rain_a_psd::Float64 = 1.0,
    rain_lw_temp::Float64 = 250.0,
    radrai_lwr::Float64 = 5.0,
    radrai_upr::Float64 = 5000.0,
    nsize_rai::Int = 64,
)
    @printf(
        "1. Computing snow band tables (Yang et al. 2013, habit=%s, Rough%03d, a=%.1f, T_LW=%.1f K)...\n",
        snow_habit,
        snow_roughness,
        snow_a_psd,
        snow_lw_temp,
    )
    sno = compute_snow_band_tables(
        data_dir;
        habit = snow_habit,
        roughness = snow_roughness,
        a_psd = snow_a_psd,
        lw_planck_temp = snow_lw_temp,
    )
    @printf(
        "   Snow r_eff range from PSD integration: [%.2f, %.2f] μm (%d points)\n",
        first(sno.r_eff),
        last(sno.r_eff),
        length(sno.r_eff),
    )

    @printf(
        "2. Computing rain band tables (Wiscombe Mie + Segelstein 1981, a=%.1f, T_LW=%.1f K)...\n",
        rain_a_psd,
        rain_lw_temp,
    )
    rai = compute_rain_band_tables(
        data_dir;
        a_psd = rain_a_psd,
        lw_planck_temp = rain_lw_temp,
    )
    @printf(
        "   Rain r_eff range from PSD integration: [%.2f, %.2f] μm (%d points)\n",
        first(rai.r_eff),
        last(rai.r_eff),
        length(rai.r_eff),
    )

    println("3. Resampling onto log-uniform effective radius grids...")
    r_sno_lw, ext_sno_lw, ssa_sno_lw, asy_sno_lw = resample_to_log_grid(
        sno.r_eff,
        sno.ext_lw,
        sno.ssa_lw,
        sno.asy_lw,
        radsno_lwr,
        radsno_upr,
        nsize_sno,
    )
    r_sno_sw, ext_sno_sw, ssa_sno_sw, asy_sno_sw = resample_to_log_grid(
        sno.r_eff,
        sno.ext_sw,
        sno.ssa_sw,
        sno.asy_sw,
        radsno_lwr,
        radsno_upr,
        nsize_sno,
    )
    r_rai_lw, ext_rai_lw, ssa_rai_lw, asy_rai_lw = resample_to_log_grid(
        rai.r_eff,
        rai.ext_lw,
        rai.ssa_lw,
        rai.asy_lw,
        radrai_lwr,
        radrai_upr,
        nsize_rai,
    )
    r_rai_sw, ext_rai_sw, ssa_rai_sw, asy_rai_sw = resample_to_log_grid(
        rai.r_eff,
        rai.ext_sw,
        rai.ssa_sw,
        rai.asy_sw,
        radrai_lwr,
        radrai_upr,
        nsize_rai,
    )

    println("4. Writing standalone precipitation NetCDF files...")
    write_precip_netcdf(
        lw_dst_nc,
        "longwave",
        RRTMGP_LW_BANDS,
        r_sno_lw,
        ext_sno_lw,
        ssa_sno_lw,
        asy_sno_lw,
        r_rai_lw,
        ext_rai_lw,
        ssa_rai_lw,
        asy_rai_lw;
        snow_habit,
        snow_roughness,
        snow_a_psd,
        snow_lw_temp,
        rain_a_psd,
        rain_lw_temp,
    )
    write_precip_netcdf(
        sw_dst_nc,
        "shortwave",
        RRTMGP_SW_BANDS,
        r_sno_sw,
        ext_sno_sw,
        ssa_sno_sw,
        asy_sno_sw,
        r_rai_sw,
        ext_rai_sw,
        ssa_rai_sw,
        asy_rai_sw;
        snow_habit,
        snow_roughness,
        snow_a_psd,
        snow_lw_temp,
        rain_a_psd,
        rain_lw_temp,
    )
    println("Done! Wrote:")
    println("  LW: $lw_dst_nc")
    println("  SW: $sw_dst_nc")
    return (; sno, rai)
end

if abspath(PROGRAM_FILE) == @__FILE__
    data_dir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "data")
    out_dir = length(ARGS) >= 2 ? ARGS[2] : @__DIR__
    snow_habit = length(ARGS) >= 3 ? ARGS[3] : "solid_column"
    snow_roughness = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 50

    suffix =
        (snow_habit == "solid_column" && snow_roughness == 50) ? "" :
        "-$(snow_habit)-$(@sprintf("Rough%03d", snow_roughness))"

    generate_rrtmgp_precip_luts(;
        data_dir,
        lw_dst_nc = joinpath(out_dir, "rrtmgp-precip-lw-bnd$(suffix).nc"),
        sw_dst_nc = joinpath(out_dir, "rrtmgp-precip-sw-bnd$(suffix).nc"),
        snow_habit,
        snow_roughness,
    )
end
