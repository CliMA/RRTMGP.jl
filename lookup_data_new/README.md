# About RRTMGP data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7988260.svg)](https://doi.org/10.5281/zenodo.7988260)

This directory contains data for use with the 
[RRTMGP k-distribution](https://github.com/earth-system-radiation/rte-rrtmgp) for computing 
absorption and scattering by gases in the the earth's atmosphere, as well an example implementation
of cloud optics and an implementation of aerosol optics following the MERRA aerosol description. 

The data have been copied and renamed from the 
[RRTMGP repository](https://github.com/earth-system-radiation/rte-rrtmgp)
at commit 74a0e09. 

Data used by the gas optics, cloud optics, and aerosol optics schemes are in the root 
directory. 

RTE+RRTMGP is distributed with some examples which are used to test implementations. 
The inputs and reference outputs for the current data are provided in `examples/`

The following files are renamed:
- `rrtmgp-gas-lw-g256.nc` => `clearsky_lw.nc`
- `rrtmgp-gas-sw-g224.nc` => `clearsky_sw.nc`
- `rrtmgp-clouds-lw.nc` => `cloudysky_lw.nc`
- `rrtmgp-clouds-sw.nc` => `cloudysky_sw.nc`
