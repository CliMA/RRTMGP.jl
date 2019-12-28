"""
    Fluxes

Compute output quantities from RTE based on spectrally-resolved flux profiles
  - This module contains an abstract class and a broadband implementation that sums over all spectral points
  - The abstract class defines the routines that extensions must implement: reduce!()
  - The intent is for users to extend it as required, using [`FluxesBroadBand`](@ref) as an example
"""
module Fluxes

using DocStringExtensions
using ..FortranIntrinsics
using ..OpticalProps

export AbstractFluxes, reduce!

"""
    AbstractFluxes{FT}

Abstract Fluxes struct
"""
abstract type AbstractFluxes{FT<:AbstractFloat} end

include("FluxesBroadBand.jl")
include("FluxesByBand.jl")

end #module
