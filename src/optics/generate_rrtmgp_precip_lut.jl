"""
    generate_rrtmgp_precip_lut.jl

Generate snow and rain optical property lookup tables for RRTMGP.jl based on:
  Feng, J., Menzel, R., & Paynter, D. (2025).
  "Toward Transparency and Consistency: An Open-Source Optics Parameterization
  for Clouds and Precipitation." Journal of Advances in Modeling Earth Systems,
  17(3), e2024MS004478. https://doi.org/10.1029/2024MS004478

This script:
1. Computes band-averaged extinction coefficient `ext` [m²/g], single-scattering
   albedo `ssa` [-], and asymmetry parameter `asy` [-] across RRTMGP's 16 longwave
   and 14 shortwave bands for:
   - Snow: Yang et al. (2013) non-spherical ice crystal database integrated over
     a Gamma/exponential PSD (default `a = 1`, matching 1M Marshall-Palmer snow).
   - Rain: Wiscombe (1980) Mie scattering with Segelstein (1981) water refractive
     indices integrated over a Gamma/exponential PSD (default `a = 1`, matching
     1M Marshall-Palmer rain).
2. Interpolates onto a log-uniform effective-radius grid:
   - Snow: `r_eff ∈ [radsno_lwr, radsno_upr]` (default `[5.0, 1100.0]` μm)
   - Rain: `r_eff ∈ [radrai_lwr, radrai_upr]` (default `[5.0, 5000.0]` μm)
3. Appends `extsno`, `ssasno`, `asysno`, `radsno_lwr`, `radsno_upr`, `radsno` and
   `extrai`, `ssarai`, `asyrai`, `radrai_lwr`, `radrai_upr`, `radrai` to the RRTMGP
   cloud lookup table NetCDF files (`rrtmgp-clouds-lw-bnd.nc` and `rrtmgp-clouds-sw-bnd.nc`).
"""

using Printf

const CM_TO_UM = 1.0e4
const M_TO_MICRON = 1.0e6
const RHO_ICE_G_M3 = 917.0e3   # g / m³ (Feng et al. 2025)
const RHO_LIQ_G_M3 = 996.0e3   # g / m³ (Feng et al. 2025)

# RRTMGP 16 Longwave bands [cm⁻¹]
const RRTMGP_LW_BANDS = [
    10.0 250.0;
    250.0 500.0;
    500.0 630.0;
    630.0 700.0;
    700.0 820.0;
    820.0 980.0;
    980.0 1080.0;
    1080.0 1180.0;
    1180.0 1390.0;
    1390.0 1480.0;
    1480.0 1800.0;
    1800.0 2080.0;
    2080.0 2250.0;
    2250.0 2390.0;
    2390.0 2680.0;
    2680.0 3250.0
]

# RRTMGP 14 Shortwave bands [cm⁻¹]
const RRTMGP_SW_BANDS = [
    820.0 2680.0;
    2680.0 3250.0;
    3250.0 4000.0;
    4000.0 4650.0;
    4650.0 5150.0;
    5150.0 6150.0;
    6150.0 7700.0;
    7700.0 8050.0;
    8050.0 12850.0;
    12850.0 16000.0;
    16000.0 22650.0;
    22650.0 29000.0;
    29000.0 38000.0;
    38000.0 50000.0
]

# ---------------------------------------------------------------------------
# 1. 1D Linear Interpolation with Linear Extrapolation (matching scipy interp1d)
# ---------------------------------------------------------------------------

"""
    interp1d_scalar(x::AbstractVector{Float64}, y::AbstractVector{Float64}, xq::Float64)

1D linear interpolation with linear extrapolation at boundaries, assuming `x` is
strictly increasing.
"""
@inline function interp1d_scalar(
    x::AbstractVector{Float64},
    y::AbstractVector{Float64},
    xq::Float64,
)
    n = length(x)
    @assert n >= 2
    if xq <= x[1]
        idx = 1
    elseif xq >= x[n]
        idx = n - 1
    else
        idx = clamp(searchsortedlast(x, xq), 1, n - 1)
    end
    @inbounds begin
        x0, x1 = x[idx], x[idx + 1]
        t = (xq - x0) / (x1 - x0)
        return muladd(t, y[idx + 1] - y[idx], y[idx])
    end
end

function interp1d(
    x::AbstractVector{Float64},
    y::AbstractVector{Float64},
    xq::AbstractVector{Float64},
)
    return [interp1d_scalar(x, y, xqi) for xqi in xq]
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
interpolate to `wavelength_um` [μm].
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
    n_interp = interp1d(freq, n_re, wavelength_um)
    k_interp = interp1d(freq, k_im, wavelength_um)
    return ComplexF64.(n_interp, .-k_interp)
end

"""
    create_size_list(var_min::Float64, var_max::Float64, n::Real = 10)

Quasi-logarithmic size grid generator from Feng et al. (2025) (`spec_util.create_list`).
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
        # Note: in Python, D[n] is 0-indexed (the (n+1)-th element of D_buf)
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
    freq = Float64[]
    d = Float64[]
    v = Float64[]
    s = Float64[]
    ext = Float64[]
    ssa = Float64[]
    asy = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            parts = split(line)
            length(parts) < 7 && continue
            push!(freq, parse(Float64, parts[1]))
            push!(d, parse(Float64, parts[2]))
            push!(v, parse(Float64, parts[3]))
            push!(s, parse(Float64, parts[4]))
            push!(ext, parse(Float64, parts[5]))
            push!(ssa, parse(Float64, parts[6]))
            push!(asy, parse(Float64, parts[7]))
        end
    end
    return freq, d, s, v, ext, ssa, asy
end

"""
    read_yang_ice_library(data_dir::AbstractString; habit="solid_column", roughness=50)

Read and merge the MIR and FIR `isca.dat` tables from Yang et al. (2013).
Returns `(wavenum_table, r_table, d_table, s_table, v_table, ext_cs, sca_cs, asy_table)`
where `ext_cs`, `sca_cs`, and `asy_table` have shape `(nwav, nsize)`.
"""
function read_yang_ice_library(
    data_dir::AbstractString;
    habit::AbstractString = "solid_column",
    roughness::Int = 50,
)
    rgh_str = @sprintf("Rough%03d", roughness)
    mir_path = joinpath(data_dir, "MIR", habit, rgh_str, "isca.dat")
    fir_path = joinpath(data_dir, "FIR", habit, rgh_str, "isca.dat")

    freq_mir, d_mir, s_mir, v_mir, ext_mir, ssa_mir, asy_mir =
        read_isca(mir_path)
    freq_fir, d_fir, s_fir, v_fir, ext_fir, ssa_fir, asy_fir =
        read_isca(fir_path)

    freq = vcat(freq_mir, freq_fir)
    d = vcat(d_mir, d_fir)
    s = vcat(s_mir, s_fir)
    v = vcat(v_mir, v_fir)
    ext_eff = vcat(ext_mir, ext_fir)
    ssa_all = vcat(ssa_mir, ssa_fir)
    asy_all = vcat(asy_mir, asy_fir)

    ext_cs_all = ext_eff .* s
    sca_cs_all = ssa_all .* ext_cs_all

    # Unique sorted frequencies and sizes (matching numpy.unique(..., return_index=True))
    freq_table = sort(unique(freq))
    d_table = sort(unique(d))
    nfreq = length(freq_table)
    nsize = length(d_table)

    s_table = zeros(Float64, nsize)
    v_table = zeros(Float64, nsize)
    for (j, dj) in enumerate(d_table)
        idx = findfirst(==(dj), d)
        s_table[j] = s[idx]
        v_table[j] = v[idx]
    end

    # Reshape (nsize, nfreq) in column-major (since size varies fastest in each wavelength block)
    ext_mat = permutedims(reshape(ext_cs_all, nsize, nfreq)) # (nfreq, nsize)
    sca_mat = permutedims(reshape(sca_cs_all, nsize, nfreq))
    asy_mat = permutedims(reshape(asy_all, nsize, nfreq))

    wavenum_unsorted = CM_TO_UM ./ freq_table
    iwav = sortperm(wavenum_unsorted)
    wavenum_table = wavenum_unsorted[iwav]

    ext_cs_table = ext_mat[iwav, :]
    sca_cs_table = sca_mat[iwav, :]
    asy_table = asy_mat[iwav, :]
    r_table = 0.75 .* v_table ./ d_table

    return wavenum_table,
    r_table,
    d_table,
    s_table,
    v_table,
    ext_cs_table,
    sca_cs_table,
    asy_table
end

# ---------------------------------------------------------------------------
# 5. Gamma PSD Integration and Spectral Band Averaging (Feng et al. 2025)
# ---------------------------------------------------------------------------

"""
    gamma_psd_integrate(a, d_char, dr, d_hres, v_hres, s_hres, ext_hres, sca_hres, asy_hres, rau)

Integrate single-particle optical properties over a Gamma PSD:
  `f(D) ∝ D^(a - 1) * exp(-a * D / d_char)`
Returns `(r_eff, s_mean, v_mean, ext, sca, ssa, asy)`.
"""
function gamma_psd_integrate(
    a::Float64,
    d_char::Float64,
    dr::Float64,
    d_hres::Vector{Float64},
    v_hres::Vector{Float64},
    s_hres::Vector{Float64},
    ext_hres::Matrix{Float64}, # (nwav, nr)
    sca_hres::Matrix{Float64}, # (nwav, nr)
    asy_hres::Matrix{Float64}, # (nwav, nr)
    rau::Float64,
)
    b = d_char / a
    nr = length(d_hres)
    nwav = size(ext_hres, 1)

    # Unnormalized gamma weights w_j = (d_j / b)^(a - 1) * exp(-d_j / b) * dr
    # Any constant factor (1 / (b * Gamma(a))) cancels out in all ratios below.
    w = Vector{Float64}(undef, nr)
    int_f = 0.0
    int_fv = 0.0
    int_fs = 0.0
    @inbounds for j in 1:nr
        dj = d_hres[j]
        wj = (dj / b)^(a - 1.0) * exp(-dj / b) * dr
        w[j] = wj
        int_f += wj
        int_fv += wj * v_hres[j]
        int_fs += wj * s_hres[j]
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
            sc = max(sca_hres[iw, j], 0.0)
            as = max(asy_hres[iw, j], 0.0)
            sum_ext += wj * max(ext_hres[iw, j], 0.0)
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
        src_sum = sum(@view source[idx])
        @inbounds for ir in 1:nr
            e_sum = 0.0
            s_sum = 0.0
            g_sum = 0.0
            for iw in idx
                sw = source[iw]
                e_iw = ext[ir, iw]
                sc_iw = e_iw * ssa[ir, iw]
                e_sum += sw * e_iw
                s_sum += sw * sc_iw
                g_sum += sw * asy[ir, iw]
            end
            ext_val = max(e_sum / src_sum, 0.0)
            sca_val = clamp(s_sum / src_sum, 0.0, ext_val)
            ssa_val = ext_val > 0.0 ? clamp(sca_val / ext_val, 0.0, 1.0) : 0.0
            asy_val = clamp(g_sum / src_sum, 0.0, 1.0)
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
    compute_snow_band_tables(data_dir; habit="solid_column", roughness=50, a_psd=1.0)

Compute LW and SW band tables `(r_eff, ext_lw, ssa_lw, asy_lw, ext_sw, ssa_sw, asy_sw)`
for snow/ice particles using Yang et al. (2013) and Feng et al. (2025).
"""
function compute_snow_band_tables(
    data_dir::AbstractString;
    habit::AbstractString = "solid_column",
    roughness::Int = 50,
    a_psd::Float64 = 1.0,
)
    wavenum_in, _, d_in, s_in, v_in, ext_cs_in, sca_cs_in, asy_in =
        read_yang_ice_library(data_dir; habit, roughness)

    r_grid = create_size_list(1.0, 10000.0, 20)
    nr = length(r_grid)
    nwav_in = length(wavenum_in)

    r_out = zeros(Float64, nr)
    ext_inres = zeros(Float64, nr, nwav_in)
    ssa_inres = zeros(Float64, nr, nwav_in)
    asy_inres = zeros(Float64, nr, nwav_in)

    cbrt_v_in = cbrt.(v_in)
    sqrt_s_in = sqrt.(s_in)
    asy_sca_in = asy_in .* sca_cs_in

    for i in 1:nr
        d_char = 2.0 * r_grid[i]
        dr = min(d_char / 100.0, 1.0)
        d_hres = collect(dr:dr:(d_char * 3.5 - 1.0e-12))
        nh = length(d_hres)

        v_hres = interp1d(d_in, cbrt_v_in, d_hres) .^ 3
        s_hres = interp1d(d_in, sqrt_s_in, d_hres) .^ 2

        ext_hres = zeros(Float64, nwav_in, nh)
        sca_hres = zeros(Float64, nwav_in, nh)
        asy_hres = zeros(Float64, nwav_in, nh)
        for iw in 1:nwav_in
            ev = @view ext_cs_in[iw, :]
            sv = @view sca_cs_in[iw, :]
            asv = @view asy_sca_in[iw, :]
            @inbounds for j in 1:nh
                sh = s_hres[j]
                e_val = interp1d_scalar(s_in, ev, sh)
                s_val = max(interp1d_scalar(s_in, sv, sh), 0.0)
                as_val =
                    s_val > 0.0 ?
                    max(interp1d_scalar(s_in, asv, sh) / s_val, 0.0) : 0.0
                ext_hres[iw, j] = e_val
                sca_hres[iw, j] = s_val
                asy_hres[iw, j] = as_val
            end
        end

        r_eff, _, _, ext_i, _, ssa_i, asy_i = gamma_psd_integrate(
            a_psd,
            d_char,
            dr,
            d_hres,
            v_hres,
            s_hres,
            ext_hres,
            sca_hres,
            asy_hres,
            RHO_ICE_G_M3,
        )
        r_out[i] = r_eff
        ext_inres[i, :] .= ext_i
        ssa_inres[i, :] .= ssa_i
        asy_inres[i, :] .= asy_i
    end

    # Interpolate to LW wavenumber grid (10:1:3249 cm⁻¹) and Planck-average
    wn_lw = collect(10.0:1.0:3249.0)
    src_lw = planck_wn(wn_lw, 250.0)
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

    # Interpolate to SW wavenumber grid (820:10:49990 cm⁻¹) and Edwards-Slingo thick-average
    wn_sw = collect(820.0:10.0:49990.0)
    wn_sol, sol_flux = read_solar_spectrum(joinpath(data_dir, "solar_flux.csv"))
    src_sw = interp1d(wn_sol, sol_flux, wn_sw)
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
            e_raw = interp1d_scalar(wn_in, e_in, wn)
            es_raw = interp1d_scalar(wn_in, es_in, wn)
            esg_raw = interp1d_scalar(wn_in, esg_in, wn)
            s_raw = e_raw != 0.0 ? es_raw / e_raw : 0.0
            g_raw = es_raw != 0.0 ? esg_raw / es_raw : 0.0
            ext_out[i, j] = max(e_raw, 0.0)
            ssa_out[i, j] = clamp(s_raw, 0.0, 1.0)
            asy_out[i, j] = clamp(g_raw, 0.0, 1.0)
        end
    end
    return ext_out, ssa_out, asy_out
end

"""
    compute_rain_band_tables(data_dir; a_psd=1.0, r_max_char=6000.0, dwn_lw=10.0, dwn_sw=50.0)

Compute LW and SW band tables `(r_eff, ext_lw, ssa_lw, asy_lw, ext_sw, ssa_sw, asy_sw)`
for liquid/rain drops using exact Wiscombe (1980) Mie scattering and Feng et al. (2025).
"""
function compute_rain_band_tables(
    data_dir::AbstractString;
    a_psd::Float64 = 1.0,
    r_max_char::Float64 = 6000.0,
    dwn_lw::Float64 = 5.0,
    dwn_sw::Float64 = 50.0,
)
    r_grid = create_size_list(0.5, r_max_char, 5)

    wn_lw = collect(10.0:dwn_lw:3249.0)
    src_lw = planck_wn(wn_lw, 250.0)
    m_lw = read_refractive_index(
        CM_TO_UM ./ wn_lw,
        joinpath(data_dir, "refractive_index"),
    )

    wn_sw = collect(820.0:dwn_sw:49990.0)
    wn_sol, sol_flux = read_solar_spectrum(joinpath(data_dir, "solar_flux.csv"))
    src_sw = interp1d(wn_sol, sol_flux, wn_sw)
    m_sw = read_refractive_index(
        CM_TO_UM ./ wn_sw,
        joinpath(data_dir, "refractive_index"),
    )

    r_out, ext_lw_hres, ssa_lw_hres, asy_lw_hres =
        _compute_mie_psd_grid(r_grid, a_psd, wn_lw, m_lw)
    _, ext_sw_hres, ssa_sw_hres, asy_sw_hres =
        _compute_mie_psd_grid(r_grid, a_psd, wn_sw, m_sw)

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

    return (; r_eff = r_out, ext_lw, ssa_lw, asy_lw, ext_sw, ssa_sw, asy_sw)
end

function _compute_mie_psd_grid(
    r_grid::Vector{Float64},
    a_psd::Float64,
    wavenum::Vector{Float64},
    m_vec::Vector{ComplexF64},
)
    # 1. Precompute exact single-particle Mie efficiencies (qext, qsca, g) on a
    #    1D log-spaced droplet radius table covering [0.05 μm, 3.5 * maximum(r_grid)],
    #    then integrate over each Gamma PSD (avoiding redundant Mie evaluations).
    r_max_single = 3.5 * maximum(r_grid)
    n_single = 180
    r_in = exp.(range(log(0.05), log(r_max_single), length = n_single))
    nwav = length(wavenum)

    qext_in = zeros(Float64, nwav, n_single)
    qsca_in = zeros(Float64, nwav, n_single)
    asy_in = zeros(Float64, nwav, n_single)

    nmed = 1.0003
    Threads.@threads for iw in 1:nwav
        D_buf = Vector{ComplexF64}(undef, 256)
        m_w = m_vec[iw]
        wn = wavenum[iw]
        @inbounds for j in 1:n_single
            x = 2.0 * π * r_in[j] * wn / CM_TO_UM * nmed
            qe, qs, g = mie_scalar(m_w, x, D_buf)
            qext_in[iw, j] = max(qe, 0.0)
            qsca_in[iw, j] = max(qs, 0.0)
            asy_in[iw, j] = max(g, 0.0)
        end
    end

    # 2. Integrate over each Gamma PSD
    nr = length(r_grid)
    r_out = zeros(Float64, nr)
    ext_out = zeros(Float64, nr, nwav)
    ssa_out = zeros(Float64, nr, nwav)
    asy_out = zeros(Float64, nr, nwav)

    log_r_in = log.(r_in)
    qasy_in = asy_in .* qsca_in

    for i in 1:nr
        d_char = 2.0 * r_grid[i]
        dr = d_char / 10.0
        d_hres = collect(dr:dr:(d_char * 3.5 - 1.0e-12))
        r_hres = 0.5 .* d_hres
        log_r_hres = log.(r_hres)
        nh = length(r_hres)

        v_hres = @. (4.0 / 3.0) * π * r_hres^3
        s_hres = @. π * r_hres^2

        ext_hres = zeros(Float64, nwav, nh)
        sca_hres = zeros(Float64, nwav, nh)
        asy_hres = zeros(Float64, nwav, nh)

        for iw in 1:nwav
            qe_row = @view qext_in[iw, :]
            qs_row = @view qsca_in[iw, :]
            qg_row = @view qasy_in[iw, :]
            @inbounds for j in 1:nh
                lr = log_r_hres[j]
                sh = s_hres[j]
                qe_val = max(interp1d_scalar(log_r_in, qe_row, lr), 0.0)
                qs_val = max(interp1d_scalar(log_r_in, qs_row, lr), 0.0)
                qg_val = max(interp1d_scalar(log_r_in, qg_row, lr), 0.0)
                ext_hres[iw, j] = qe_val * sh
                sca_hres[iw, j] = qs_val * sh
                asy_hres[iw, j] = qs_val > 0.0 ? qg_val / qs_val : 0.0
            end
        end

        r_eff, _, _, ext_i, _, ssa_i, asy_i = gamma_psd_integrate(
            a_psd,
            d_char,
            dr,
            d_hres,
            v_hres,
            s_hres,
            ext_hres,
            sca_hres,
            asy_hres,
            RHO_LIQ_G_M3,
        )
        r_out[i] = r_eff
        ext_out[i, :] .= ext_i
        ssa_out[i, :] .= ssa_i
        asy_out[i, :] .= asy_i
    end

    return r_out, ext_out, ssa_out, asy_out
end

# ---------------------------------------------------------------------------
# 7. Resample onto Log-Uniform Effective Radius Grids for RRTMGP.jl
# ---------------------------------------------------------------------------

"""
    resample_to_log_grid(r_src, ext_src, ssa_src, asy_src, r_lwr, r_upr, nsize)

Resample band tables `(nr_src, nband)` onto a log-uniform effective radius grid:
  `r_grid[k] = r_lwr * (r_upr / r_lwr)^((k - 1) / (nsize - 1))` for `k = 1:nsize`.
Interpolates in `log(r_eff)` so that `compute_lookup_cld_*_props` in RRTMGP.jl
can do O(1) direct index lookup in `log(r_eff)`.
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
    log_r_src = log.(r_src)
    log_r_grid = range(log(r_lwr), log(r_upr), length = nsize)
    r_grid = exp.(collect(log_r_grid))

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
                Float32(max(interp1d_scalar(log_r_src, e_col, lr), 0.0))
            ssa_lut[k, ib] =
                Float32(clamp(interp1d_scalar(log_r_src, s_col, lr), 0.0, 1.0))
            asy_lut[k, ib] =
                Float32(clamp(interp1d_scalar(log_r_src, g_col, lr), 0.0, 1.0))
        end
    end

    return Float32.(r_grid), ext_lut, ssa_lut, asy_lut
end

# ---------------------------------------------------------------------------
# 8. Append Snow & Rain Variables to RRTMGP Cloud NetCDF Files
# ---------------------------------------------------------------------------

"""
    append_precip_to_rrtmgp_netcdf(
        src_nc_path, dst_nc_path, λ_str,
        r_sno, ext_sno, ssa_sno, asy_sno,
        r_rai, ext_rai, ssa_rai, asy_rai,
    )

Copy `src_nc_path` to `dst_nc_path` and append snow (`extsno`, `ssasno`, `asysno`,
`radsno_lwr`, `radsno_upr`, `radsno`) and rain (`extrai`, `ssarai`, `asyrai`,
`radrai_lwr`, `radrai_upr`, `radrai`) lookup tables.
Uses `NCDatasets` when available in the active Julia environment, or falls back
to Python `netCDF4`.
"""
function append_precip_to_rrtmgp_netcdf(
    src_nc_path::AbstractString,
    dst_nc_path::AbstractString,
    λ_str::AbstractString,
    r_sno::Vector{Float32},
    ext_sno::Matrix{Float32}, # (nsize_sno, nband) in Julia column-major
    ssa_sno::Matrix{Float32},
    asy_sno::Matrix{Float32},
    r_rai::Vector{Float32},
    ext_rai::Matrix{Float32}, # (nsize_rai, nband) in Julia column-major
    ssa_rai::Matrix{Float32},
    asy_rai::Matrix{Float32},
)
    cp(src_nc_path, dst_nc_path; force = true)
    if Base.find_package("NCDatasets") !== nothing
        @eval import NCDatasets
        Base.invokelatest(
            _append_ncdatasets!,
            dst_nc_path,
            λ_str,
            r_sno,
            ext_sno,
            ssa_sno,
            asy_sno,
            r_rai,
            ext_rai,
            ssa_rai,
            asy_rai,
        )
    else
        _append_python_netcdf4!(
            dst_nc_path,
            λ_str,
            r_sno,
            ext_sno,
            ssa_sno,
            asy_sno,
            r_rai,
            ext_rai,
            ssa_rai,
            asy_rai,
        )
    end
    return dst_nc_path
end

function _append_ncdatasets!(
    dst_nc_path,
    λ_str,
    r_sno,
    ext_sno,
    ssa_sno,
    asy_sno,
    r_rai,
    ext_rai,
    ssa_rai,
    asy_rai,
)
    NCD = @eval NCDatasets
    NCD.NCDataset(dst_nc_path, "a") do ds
        NCD.defDim(ds, "nsize_sno", length(r_sno))
        NCD.defDim(ds, "nsize_rai", length(r_rai))

        NCD.defVar(
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
        NCD.defVar(
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
        NCD.defVar(
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
        NCD.defVar(
            ds,
            "extsno",
            ext_sno,
            ("nsize_sno", "nband"),
            attrib = Dict(
                "standard_name" => "extsno",
                "description" => "Band LUT $(λ_str) snow extinction coefficient (Feng et al. 2025 / Yang et al. 2013)",
                "units" => "m2/g",
            ),
        )
        NCD.defVar(
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
        NCD.defVar(
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

        NCD.defVar(
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
        NCD.defVar(
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
        NCD.defVar(
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
        NCD.defVar(
            ds,
            "extrai",
            ext_rai,
            ("nsize_rai", "nband"),
            attrib = Dict(
                "standard_name" => "extrai",
                "description" => "Band LUT $(λ_str) rain extinction coefficient (Feng et al. 2025 Mie)",
                "units" => "m2/g",
            ),
        )
        NCD.defVar(
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
        NCD.defVar(
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
end

function _append_python_netcdf4!(
    dst_nc_path,
    λ_str,
    r_sno,
    ext_sno,
    ssa_sno,
    asy_sno,
    r_rai,
    ext_rai,
    ssa_rai,
    asy_rai,
)
    # Write binary scratch file and invoke python3 netCDF4 to append in C-order (nband, nsize)
    tmp_bin = dst_nc_path * ".tmp.bin"
    open(tmp_bin, "w") do io
        write(
            io,
            Int32(length(r_sno)),
            Int32(length(r_rai)),
            Int32(size(ext_sno, 2)),
        )
        write(io, r_sno)
        # Note: in Python netCDF4 (C-order), ("nband", "nsize_sno") corresponds to
        # Julia column-major ("nsize_sno", "nband")!
        write(io, ext_sno)
        write(io, ssa_sno)
        write(io, asy_sno)
        write(io, r_rai)
        write(io, ext_rai)
        write(io, ssa_rai)
        write(io, asy_rai)
    end
    py_code = """
import netCDF4 as nc, numpy as np, sys
nc_path, bin_path, lam_str = sys.argv[1], sys.argv[2], sys.argv[3]
with open(bin_path, "rb") as f:
    nsize_sno, nsize_rai, nband = np.fromfile(f, dtype=np.int32, count=3)
    r_sno = np.fromfile(f, dtype=np.float32, count=nsize_sno)
    ext_sno = np.fromfile(f, dtype=np.float32, count=nsize_sno*nband).reshape((nband, nsize_sno))
    ssa_sno = np.fromfile(f, dtype=np.float32, count=nsize_sno*nband).reshape((nband, nsize_sno))
    asy_sno = np.fromfile(f, dtype=np.float32, count=nsize_sno*nband).reshape((nband, nsize_sno))
    r_rai = np.fromfile(f, dtype=np.float32, count=nsize_rai)
    ext_rai = np.fromfile(f, dtype=np.float32, count=nsize_rai*nband).reshape((nband, nsize_rai))
    ssa_rai = np.fromfile(f, dtype=np.float32, count=nsize_rai*nband).reshape((nband, nsize_rai))
    asy_rai = np.fromfile(f, dtype=np.float32, count=nsize_rai*nband).reshape((nband, nsize_rai))

with nc.Dataset(nc_path, "a") as ds:
    ds.createDimension("nsize_sno", int(nsize_sno))
    ds.createDimension("nsize_rai", int(nsize_rai))
    for name, val, desc in [
        ("radsno_lwr", r_sno[0], "Snow particle effective radius lower bound for band LUT log-interpolation"),
        ("radsno_upr", r_sno[-1], "Snow particle effective radius upper bound for band LUT log-interpolation"),
        ("radrai_lwr", r_rai[0], "Rain particle effective radius lower bound for band LUT log-interpolation"),
        ("radrai_upr", r_rai[-1], "Rain particle effective radius upper bound for band LUT log-interpolation"),
    ]:
        v = ds.createVariable(name, "f4", ())
        v.standard_name = name
        v.description = desc
        v.units = "microns"
        v[:] = val

    vs = ds.createVariable("radsno", "f4", ("nsize_sno",))
    vs.standard_name = "radsno"
    vs.description = "Snow particle effective radius log-uniform grid points"
    vs.units = "microns"
    vs[:] = r_sno

    vr = ds.createVariable("radrai", "f4", ("nsize_rai",))
    vr.standard_name = "radrai"
    vr.description = "Rain particle effective radius log-uniform grid points"
    vr.units = "microns"
    vr[:] = r_rai

    for name, arr, units, prop in [
        ("extsno", ext_sno, "m2/g", "snow extinction coefficient"),
        ("ssasno", ssa_sno, "unitless", "snow single scattering albedo"),
        ("asysno", asy_sno, "unitless", "snow asymmetry parameter"),
    ]:
        v = ds.createVariable(name, "f4", ("nband", "nsize_sno"))
        v.standard_name = name
        v.description = f"Band LUT {lam_str} {prop} (Feng et al. 2025 / Yang et al. 2013)"
        v.units = units
        v[:, :] = arr

    for name, arr, units, prop in [
        ("extrai", ext_rai, "m2/g", "rain extinction coefficient"),
        ("ssarai", ssa_rai, "unitless", "rain single scattering albedo"),
        ("asyrai", asy_rai, "unitless", "rain asymmetry parameter"),
    ]:
        v = ds.createVariable(name, "f4", ("nband", "nsize_rai"))
        v.standard_name = name
        v.description = f"Band LUT {lam_str} {prop} (Feng et al. 2025 Mie)"
        v.units = units
        v[:, :] = arr
"""
    run(`python3 -c $py_code $dst_nc_path $tmp_bin $λ_str`)
    rm(tmp_bin; force = true)
end

"""
    generate_rrtmgp_precip_luts(;
        data_dir,
        lw_src_nc,
        sw_src_nc,
        lw_dst_nc,
        sw_dst_nc,
        snow_habit = "solid_column",
        snow_roughness = 50,
        snow_a_psd = 1.0,
        radsno_lwr = 5.0,
        radsno_upr = 1100.0,
        nsize_sno = 64,
        rain_a_psd = 1.0,
        radrai_lwr = 5.0,
        radrai_upr = 5000.0,
        nsize_rai = 64,
    )

End-to-end generator for RRTMGP longwave and shortwave cloud+precipitation lookup tables.
"""
function generate_rrtmgp_precip_luts(;
    data_dir::AbstractString,
    lw_src_nc::AbstractString,
    sw_src_nc::AbstractString,
    lw_dst_nc::AbstractString,
    sw_dst_nc::AbstractString,
    snow_habit::AbstractString = "solid_column",
    snow_roughness::Int = 50,
    snow_a_psd::Float64 = 1.0,
    radsno_lwr::Float64 = 5.0,
    radsno_upr::Float64 = 1100.0,
    nsize_sno::Int = 64,
    rain_a_psd::Float64 = 1.0,
    radrai_lwr::Float64 = 5.0,
    radrai_upr::Float64 = 5000.0,
    nsize_rai::Int = 64,
)
    println(
        "1. Computing snow band tables (Yang et al. 2013, habit=$snow_habit, Rough$(snow_roughness), a=$snow_a_psd)...",
    )
    sno = compute_snow_band_tables(
        data_dir;
        habit = snow_habit,
        roughness = snow_roughness,
        a_psd = snow_a_psd,
    )
    @printf(
        "   Snow r_eff range from PSD integration: [%.2f, %.2f] μm (%d points)\n",
        first(sno.r_eff),
        last(sno.r_eff),
        length(sno.r_eff)
    )

    println(
        "2. Computing rain band tables (Wiscombe Mie + Segelstein 1981, a=$rain_a_psd)...",
    )
    rai = compute_rain_band_tables(data_dir; a_psd = rain_a_psd)
    @printf(
        "   Rain r_eff range from PSD integration: [%.2f, %.2f] μm (%d points)\n",
        first(rai.r_eff),
        last(rai.r_eff),
        length(rai.r_eff)
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

    println("4. Appending snow and rain tables to NetCDF files...")
    append_precip_to_rrtmgp_netcdf(
        lw_src_nc,
        lw_dst_nc,
        "longwave",
        r_sno_lw,
        ext_sno_lw,
        ssa_sno_lw,
        asy_sno_lw,
        r_rai_lw,
        ext_rai_lw,
        ssa_rai_lw,
        asy_rai_lw,
    )
    append_precip_to_rrtmgp_netcdf(
        sw_src_nc,
        sw_dst_nc,
        "shortwave",
        r_sno_sw,
        ext_sno_sw,
        ssa_sno_sw,
        asy_sno_sw,
        r_rai_sw,
        ext_rai_sw,
        ssa_rai_sw,
        asy_rai_sw,
    )
    println("Done! Wrote:")
    println("  LW: $lw_dst_nc")
    println("  SW: $sw_dst_nc")
    return (; sno, rai)
end

if abspath(PROGRAM_FILE) == @__FILE__
    base_dir = @__DIR__
    generate_rrtmgp_precip_luts(;
        data_dir = joinpath(base_dir, "data"),
        lw_src_nc = joinpath(base_dir, "rrtmgp-clouds-lw-bnd.nc"),
        sw_src_nc = joinpath(base_dir, "rrtmgp-clouds-sw-bnd.nc"),
        lw_dst_nc = joinpath(base_dir, "rrtmgp-clouds-lw-bnd-with-precip.nc"),
        sw_dst_nc = joinpath(base_dir, "rrtmgp-clouds-sw-bnd-with-precip.nc"),
    )
end
