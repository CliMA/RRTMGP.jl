# Optics

```@meta
CurrentModule = RRTMGP.Optics
```

```@docs
AbstractOpticalProps
OneScalar
TwoStream
compute_col_gas!
compute_optical_props!
```

```@docs
compute_col_gas_kernel!
compute_interp_fractions
compute_interp_frac_temp
compute_interp_frac_press
compute_interp_frac_η
compute_τ_ssa_lw_src!
compute_τ_minor
compute_τ_rayleigh
compute_lw_planck_src!
```

```@docs
build_cloud_mask!
add_cloud_optics_2stream
compute_cld_props
pade_eval
```

```@docs
compute_optical_props_kernel!
compute_optical_props_kernel_lw!
compute_sources_gray_kernel!
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
