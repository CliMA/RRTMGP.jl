import ..AtmosphericStates: getview_p_lay, getview_t_lay, getview_rel_hum
import ..Vmrs
# TODO: finish getter methods
#! format: off

include("getters.jl")
optical_thickness_parameter(s::RRTMGPSolver) = nothing                                          # not yet used

isothermal_boundary_layer(s::RRTMGPSolver) = s.grid_params.isothermal_boundary_layer

"""
    domain_view(gp::RRTMGPGridParams, data::RRTMGPData{NVCOrder})

Return a view into the domain portion of the data, (excludes the isothermal
boundary layer if it exists)
"""
function domain_view end

@inline domain_view(_, data::RRTMGPData) = data
@inline domain_view(gp::RRTMGPGridParams, data::RRTMGPData{NVCOrder}) =
    domain_view(gp.isothermal_boundary_layer, data)

@inline domain_view(s::RRTMGPSolver, data::RRTMGPData{NVCOrder}) =
    domain_view(isothermal_boundary_layer(s), data)

function domain_view(isothermal_boundary_layer::Bool, data::RRTMGPData{NVCOrder})
    nlay = size(parent(data), 2)
    if isothermal_boundary_layer
        RRTMGPData{NVCOrder}(view(parent(data), :, 1:nlay-1, :))
    else
        RRTMGPData{NVCOrder}(view(parent(data), :, 1:nlay, :))
    end
end

function domain_view(isothermal_boundary_layer::Bool, data::RRTMGPData{VCOrder})
    nlay = size(parent(data), 1)
    if isothermal_boundary_layer
        RRTMGPData{VCOrder}(view(parent(data), 1:nlay-1, :))
    else
        RRTMGPData{VCOrder}(view(parent(data), 1:nlay, :))
    end
end

"""
    aerosol_names()

Return a vector containing aerosol names used in the `AerosolState`.
"""
aerosol_names() = [
    "dust1",
    "ss1",
    "so4",
    "bcpi",
    "bcpo",
    "ocpi",
    "ocpo",
    "dust2",
    "dust3",
    "dust4",
    "dust5",
    "ss2",
    "ss3",
    "ss4",
    "ss5",
]

for (i, aname) in enumerate(aerosol_names())
    sname = Symbol(aname)
    (occursin("dust", aname) || occursin("ss", aname)) || continue

    # Double interpolation needed here: 1 for @eval, second for doc string
    @eval begin
        """
            radius_$($(aname))(s::RRTMGPSolver)

        Returns the aerosol radius for species $($(aname)).
        """
        $(Symbol(:radius_, sname))(s::RRTMGPSolver) = view(s.as.aerosol_state.aero_size, $i, :, :) # center_$(name)_radius

        """
            column_mass_density_$($(aname))(s::RRTMGPSolver)

        Returns the column aerosol mass density (volume mixing ratio)
        for species $($(aname)).
        """
        $(Symbol(:column_mass_density_, sname))(s::RRTMGPSolver) = view(s.as.aerosol_state.aero_mass, $i, :, :) # center_$(name)_radius
    end
end

"""
    gas_names_sw()

Return a vector containing the gas names in the shortwave lookup tables.

This should be the same list of gases returned from the following code,
and is tested in `lookup_tables`
```julia
function gas_names_sw_from_artifacts()
   artifact(t, b, n) =
       NC.Dataset(RRTMGP.ArtifactPaths.get_lookup_filename(t, b)) do ds
           getproperty(RRTMGP.LookUpTables, n)(ds, Float64, Array)
       end
   _, idx_gases_sw = artifact(:gas, :sw, :LookUpSW)
   return keys(idx_gases_sw)
end
```
"""
gas_names_sw() = [
    "h2o",
    "cfc11",
    "h2o_self",
    "co2",
    "cfc12",
    "hfc134a",
    "cfc22",
    "ch4",
    "hfc23",
    "ccl4",
    "hfc143a",
    "co",
    "no2",
    "n2",
    "o2",
    "o3",
    "h2o_frgn",
    "hfc32",
    "n2o",
    "cf4",
    "hfc125"
]

for gname in gas_names_sw()
    sname = Symbol(gname)
    # Double interpolation needed here: 1 for @eval, second for doc string
    @eval begin
        """
            volume_mixing_ratio_$($(gname))(s::RRTMGPSolver)

        Returns the volume moxing ratio for $($(gname)).
        """
        $(Symbol(:volume_mixing_ratio_, sname))(s::RRTMGPSolver) =
            return get_vmr(s.as.vmr, s.lookups.idx_gases_sw[gname])
    end
end

get_vmr(vmr::Vmrs.VmrGM, idx_gase_sw) = view(vmr.vmr, idx_gase_sw)
get_vmr(vmr::Vmrs.Vmr, idx_gases_sw) =
    view(vmr.vmr, idx_gase_sw, :, :)

# #! format: on
