# Dataset files

# This function generates the file names for lookup table files, 
# atmospheric state files and reference data file for comparing
# results for test cases.
# The following options are currently supported:
# ftype ∈ (:lookup_tables, :atmos_state, :comparison)
# optics_type ∈ (:clearsky, :cloudy)
# opc ∈ (:OneScalar, :TwoStream, nothing)
# λ ∈ (:lw, :sw, nothing)
# flux_up_dn ∈ (:flux_up, :flux_dn, nothing)

function get_ref_filename(
    ftype,
    optics_type;
    λ = nothing,
    opc = nothing,
    flux_up_dn = nothing,
)
    @assert ftype ∈ (:lookup_tables, :atmos_state, :comparison) &&
            optics_type ∈ (:clearsky, :cloudysky)

    fname = String(optics_type)

    if ftype == :lookup_tables
        @assert λ ∈ (:lw, :sw)
        fname *= "_" * String(λ)
    elseif ftype == :atmos_state
        fname *= "_as"
    else #ftype == :comparison
        @assert opc ∈ (:OneScalar, :TwoStream) &&
                λ ∈ (:lw, :sw) &&
                flux_up_dn ∈ (:flux_up, :flux_dn)
        fname *= "_" * String(λ) * "_" * String(flux_up_dn) * "_" * String(opc)
    end

    return joinpath(
        artifact"RRTMGPReferenceData",
        "RRTMGPReferenceData",
        String(ftype),
        fname * ".nc",
    )
end
