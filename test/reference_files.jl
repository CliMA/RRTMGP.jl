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

function get_ref_filename(ftype, optics_type; λ = nothing, opc = nothing, flux_up_dn = nothing, lut_pade = nothing)
    @assert ftype ∈ (:lookup_tables, :atmos_state, :comparison) && optics_type ∈ (:clearsky, :cloudysky)

    fname = String(optics_type)

    if ftype == :lookup_tables
        @assert λ ∈ (:lw, :sw)
        fname *= "_" * String(λ)
    elseif ftype == :atmos_state
        fname *= "_as"
    else #ftype == :comparison
        if λ isa Symbol
            @assert λ ∈ (:lw, :sw)
            fname *= "_" * String(λ)
        end

        if flux_up_dn isa Symbol
            @assert flux_up_dn ∈ (:flux_up, :flux_dn)
            fname *= "_" * String(flux_up_dn)
        end
        if opc isa Symbol
            @assert opc ∈ (:OneScalar, :TwoStream)
            fname *= "_" * String(opc)
        end
        if lut_pade isa Symbol
            @assert lut_pade ∈ (:lut, :pade)
            fname *= "_" * String(lut_pade)
        end
    end
    if ftype == :lookup_tables
        dir = RRTMGP.lookup_data()
    else
        dir = joinpath(artifact"RRTMGPReferenceData", "RRTMGPReferenceData", String(ftype))
    end

    return joinpath(dir, fname * ".nc")
end
