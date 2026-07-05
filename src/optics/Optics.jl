module Optics

using DocStringExtensions
using Adapt
using Random
import ClimaComms

using ..Vmrs
import ..pow_fast
using ..LookUpTables
using ..AtmosphericStates
import ..RRTMGPGridParams
using ..Sources
using ..AngularDiscretizations
import ..Parameters as RP

export AbstractOpticalProps,
    OneScalar,
    TwoStream,
    compute_col_gas!,
    compute_relative_humidity!,
    compute_optical_props!

"""
    AbstractOpticalProps

Optical properties for one scalar and two stream calculations.
"""
abstract type AbstractOpticalProps end

# Optical-property containers (OneScalar/TwoStream) and their constructors
include("optical_props.jl")
# Column dry-air amounts and relative humidity (thermodynamic helpers)
include("column_amounts.jl")
# compute_optical_props! dispatch: gray/spectral × longwave/shortwave
include("compute_optical_props.jl")

# Interpolation helpers and the per-domain optics kernels
include("optics_utils.jl")
include("gas_optics.jl")
include("cloud_optics.jl")
include("aerosol_optics.jl")
include("gray_optics_kernels.jl")

end
