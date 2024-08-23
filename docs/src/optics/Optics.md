# Optics

```@meta
CurrentModule = RRTMGP.Optics
```

```@docs
AbstractOpticalProps
OneScalar
TwoStream
compute_col_gas!
compute_relative_humidity!
compute_optical_props!
```

```@docs
compute_col_gas_kernel!
compute_relative_humidity_kernel!
compute_interp_frac_temp
compute_interp_frac_press
compute_interp_frac_η
compute_gas_optics
compute_τ_minor
compute_τ_rayleigh
```

```@docs
build_cloud_mask!
add_cloud_optics_2stream!
compute_lookup_cld_liq_props
compute_lookup_cld_ice_props
compute_pade_cld_props
pade_eval
```

```@docs
add_aerosol_optics_1scalar!
add_aerosol_optics_2stream!
compute_lookup_aerosol
compute_lookup_dust_props
compute_lookup_sea_salt_props
compute_lookup_sulfate_props
compute_lookup_black_carbon_props
compute_lookup_organic_carbon_props
locate_merra_size_bin
```

```@docs
compute_gray_optical_thickness_lw
compute_gray_optical_thickness_sw
```

```@docs
interp1d
interp2d
interp3d
increment_2stream
delta_scale
```
