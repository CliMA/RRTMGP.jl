# RRTMGP lookup data

This is repackaged data from v1.0 of [rte-rrtmgp](https://github.com/earth-system-radiation/rte-rrtmgp), 

> Robert Pincus, Benjamin R. Hillman, Matt Norman, fomics, & Chiel van Heerwaarden. (2019). RobertPincus/rte-rrtmgp: (Re-)release for 10.1029/2019MS001621 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.3403173

suitable for use as a Julia artifact with RRTMGP.jl

Specifically
- `rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc` => `clearsky_lw.nc`
- `rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc` => `clearsky_sw.nc`
- `extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc` => `cloudysky_lw.nc`
- `extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc` => `cloudysky_sw.nc`

Note: data for versions up to v1.6 are stored in the above repository; data for later versions are stored in a separate repository https://github.com/earth-system-radiation/rrtmgp-data.
