
# use mo_rte_kind,      only: wp, wl
# use mo_rte_config,    only: check_extents, check_values
# use mo_rte_util_array_validation,&
#                       only: extents_are, any_vals_outside
# use mo_optical_props, only: ty_optical_props,      &
#                             ty_optical_props_arry, &
#                             ty_optical_props_1scl, &
#                             ty_optical_props_2str, &
#                             ty_optical_props_nstr
# implicit none

# MERRA2/GOCART aerosol types
const merra_ntype = 7          # Number of MERRA aerosol types
const merra_aero_none = 0      # no aerosal
const merra_aero_dust = 1      # dust
const merra_aero_salt = 2      # Salt
const merra_aero_sulf = 3      # sulfate
const merra_aero_bcar_rh = 4   # black carbon, hydrophilic
const merra_aero_bcar = 5      # black carbon, hydrophobic
const merra_aero_ocar_rh = 6   # organic carbon, hydrophilic
const merra_aero_ocar = 7      # organic carbon, hydrophobic

# index identifiers for aerosol optical property tables
const ext = 1                  # extinction
const ssa = 2                  # single scattering albedo
const g   = 3                  # asymmetry parameter

# private
# -----------------------------------------------------------------------------------
abstract type AbstractAerosolOpticsMerra <: AbstractOpticalProps end

struct AerosolOpticsMerra{FTA1D,FTA2D,FTA3D,FTA4D} <: AbstractAerosolOpticsMerra
# type, extends(ty_optical_props), public :: ty_aerosol_optics_rrtmgp_merra
  # private
  #
  # Lookup table information
  #
  # Table upper and lower aerosol size (radius) bin limits (microns)
  merra_aero_bin_lims::FTA2D     # Dimensions (pair,nbin)
  # Table relative humidity values
  aero_rh::FTA1D
  #
  # The aerosol tables themselves.
  # extinction (m2/kg)
  # single scattering albedo (unitless)
  # asymmetry parameter (unitless)
  #
  aero_dust_tbl::FTA3D      # ext, ssa, g (nval, nbin, nbnd)
  aero_salt_tbl::FTA4D      # ext, ssa, g (nval, nrh, nbin, nbnd)
  aero_sulf_tbl::FTA3D      # ext, ssa, g (nval, nrh, nbnd)
  aero_bcar_tbl::FTA2D      # ext, ssa, g (nval, nbnd)
  aero_bcar_rh_tbl::FTA3D   # ext, ssa, g (nval, nrh, nbnd)
  aero_ocar_tbl::FTA2D      # ext, ssa, g (nval, nbnd)
  aero_ocar_rh_tbl::FTA3D   # ext, ssa, g (nval, nrh, nbnd)
  #
  # -----
# contains
  # generic,   public  :: load  => load_lut
  # procedure, public  :: finalize
  # procedure, public  :: aerosol_optics
  # # Internal procedures
  # procedure, private :: load_lut
end

# ------------------------------------------------------------------------------
#
# Routines to load data needed for aerosol optics calculations from lookup-tables.
#
# ------------------------------------------------------------------------------
function load_lut(this, band_lims_wvn,
                merra_aero_bin_lims, aero_rh,
                aero_dust_tbl, aero_salt_tbl, aero_sulf_tbl,
                aero_bcar_tbl, aero_bcar_rh_tbl,
                aero_ocar_tbl, aero_ocar_rh_tbl; check_extents=true)
                  # result(error_msg)

  # class(ty_aerosol_optics_rrtmgp_merra),   &
  #                             intent(inout) :: this
  # real(wp), dimension(:,:),   intent(in   ) :: band_lims_wvn ! spectral discretization
  # ! Lookup table interpolation constants
  # real(wp), dimension(:,:),   intent(in   ) :: merra_aero_bin_lims ! aerosol lut size bin limiits (pair,nbin)
  # real(wp), dimension(:),     intent(in   ) :: aero_rh       ! relative humidity LUT dimension values
  # ! LUT coefficients
  # ! Extinction, single-scattering albedo, and asymmetry parameter for aerosol types
  # real(wp), dimension(:,:,:),   intent(in)  :: aero_dust_tbl
  # real(wp), dimension(:,:,:,:), intent(in)  :: aero_salt_tbl
  # real(wp), dimension(:,:,:),   intent(in)  :: aero_sulf_tbl
  # real(wp), dimension(:,:),     intent(in)  :: aero_bcar_tbl
  # real(wp), dimension(:,:,:),   intent(in)  :: aero_bcar_rh_tbl
  # real(wp), dimension(:,:),     intent(in)  :: aero_ocar_tbl
  # real(wp), dimension(:,:,:),   intent(in)  :: aero_ocar_rh_tbl
  # character(len=128)    :: error_msg
  # -------
  #
  # Local variables
  #
  # integer               :: npair, nval, nrh, nbin, nband

  # error_msg = this.init(band_lims_wvn, name="RRTMGP aerosol optics")
  #
  # LUT coefficient dimensions
  #
  npair = size(merra_aero_bin_lims,dim=1)
  nval  = size(aero_salt_tbl,dim=1)
  nrh   = size(aero_salt_tbl,dim=2)
  nbin  = size(aero_salt_tbl,dim=3)
  nband = size(aero_salt_tbl,dim=4)
  #
  # Error checking
  #
  if (check_extents)
    extents_are(aero_dust_tbl, nval, nbin, nband) || return  "aerosol_optics.load_lut(): array aero_dust_tbl isn't consistently sized"
    extents_are(aero_salt_tbl, nval, nrh, nbin, nband) || return  "aerosol_optics.load_lut(): array aero_salt_tbl isn't consistently sized"
    extents_are(aero_sulf_tbl, nval, nrh, nband) || return  "aerosol_optics.load_lut(): array aero_sulf_tbl isn't consistently sized"
    extents_are(aero_bcar_rh_tbl, nval, nrh, nband) || return  "aerosol_optics.load_lut(): array aero_bcar_rh_tbl isn't consistently sized"
    extents_are(aero_bcar_tbl, nval, nband) || return  "aerosol_optics.load_lut(): array aero_bcar_tbl isn't consistently sized"
    extents_are(aero_ocar_rh_tbl, nval, nrh, nband) || return  "aerosol_optics.load_lut(): array aero_ocar_rh_tbl isn't consistently sized"
    extents_are(aero_ocar_tbl, nval, nband) || return  "aerosol_optics.load_lut(): array aero_ocar_tbl isn't consistently sized"
  end

  # Allocate LUT parameters
  this.merra_aero_bin_lims = zeros(npair,nbin)
  this.aero_rh = zeros(nrh)
  # Allocate LUT coefficients
  this.aero_dust_tbl = zeros(nval, nbin, nband)
  this.aero_salt_tbl = zeros(nrh, nval, nbin, nband)
  this.aero_sulf_tbl = zeros(nrh, nval, nband)
  this.aero_bcar_tbl = zeros(nval, nband)
  this.aero_bcar_rh_tbl = zeros(nrh, nval, nband)
  this.aero_ocar_tbl = zeros(nval, nband)
  this.aero_ocar_rh_tbl = zeros(nrh, nval, nband)

  # Copy LUT coefficients
  this.merra_aero_bin_lims .= merra_aero_bin_lims
  this.aero_rh             .= aero_rh
  this.aero_dust_tbl .= aero_dust_tbl
  this.aero_bcar_tbl .= aero_bcar_tbl
  this.aero_ocar_tbl .= aero_ocar_tbl

  this.aero_salt_tbl    .= freshape( aero_salt_tbl;    shape=(nrh, nval, nbin, nband), order=(2,1,3,4) )
  this.aero_sulf_tbl    .= freshape( aero_sulf_tbl;    shape=(nrh, nval,       nband), order=(2,1,3) )
  this.aero_bcar_rh_tbl .= freshape( aero_bcar_rh_tbl; shape=(nrh, nval,       nband), order=(2,1,3) )
  this.aero_ocar_rh_tbl .= freshape( aero_ocar_rh_tbl; shape=(nrh, nval,       nband), order=(2,1,3) )

  # !$acc enter data create(this)
  # !$acc            copyin(this.aero_dust_tbl, this.aero_salt_tbl, this.aero_sulf_tbl)
  # !$acc            copyin(this.aero_bcar_tbl, this.aero_bcar_rh_tbl)
  # !$acc            copyin(this.aero_ocar_tbl, this.aero_ocar_rh_tbl)
  # !$acc            copyin(this.merra_aero_bin_lims, this.aero_rh)
  # !$omp target enter data
  # !$omp map(to:this.aero_dust_tbl, this.aero_salt_tbl, this.aero_sulf_tbl)
  # !$omp map(to:this.aero_bcar_tbl, this.aero_bcar_rh_tbl)
  # !$omp map(to:this.aero_ocar_tbl, this.aero_ocar_rh_tbl)
  # !$omp map(to:this.merra_aero_bin_lims, this.aero_rh)

end
# ------------------------------------------------------------------------------
#
# Derive aerosol optical properties from provided aerosol input properties
#
# ------------------------------------------------------------------------------
#
# Compute single-scattering properties
#
function aerosol_optics(this, aero_type, aero_size, aero_mass, relhum,
                        optical_props; check_values=true)::String

  # class(ty_aerosol_optics_rrtmgp_merra), &
  #           intent(in  ) :: this
  # integer,  intent(in  ) :: aero_type(:,:)   ! MERRA2/GOCART aerosol type
  #                                            ! Dimensions: (ncol,nlay)
  #                                            ! 1 = merra_aero_dust    (dust)
  #                                            ! 2 = merra_aero_salt    (salt)
  #                                            ! 3 = merra_aero_sulf    (sulfate)
  #                                            ! 4 = merra_aero_bcar_rh (black carbon, hydrophilic)
  #                                            ! 5 = merra_aero_bcar    (black carbon, hydrophobic)
  #                                            ! 6 = merra_aero_ocar_rh (organic carbon, hydrophilic)
  #                                            ! 7 = merra_aero_ocar    (organic carbon, hydrophobic)
  # real(wp), intent(in  ) :: aero_size(:,:)   ! aerosol size (radius) for dust and sea-salt (microns)
  #                                            ! Dimensions: (ncol,nlay)
  # real(wp), intent(in  ) :: aero_mass(:,:)   ! aerosol mass column (kg/m2)
  #                                            ! Dimensions: (ncol,nlay)
  # real(wp), intent(in  ) :: relhum(:,:)      ! relative humidity (fraction, 0-1)
  #                                            ! Dimensions: (ncol,nlay)

  # class(ty_optical_props_arry), &
  #           intent(inout) :: optical_props
                                             # Dimensions: (ncol,nlay,nbnd)

  # character(len=128)      :: error_msg
  # ------- Local -------
  # logical(wl), dimension(size(aero_type,1), size(aero_type,2)) :: aeromsk
  # real(wp),    dimension(size(aero_type,1), size(aero_type,2), size(this.aero_dust_tbl,3)) :: &
  #              atau, ataussa, ataussag
  # integer  :: ncol, nlay, npair, nbin, nrh, nval, nbnd
  # integer  :: icol, ilay, ibnd, ibin
  # scalars for total tau, tau*ssa
  # real(wp) :: tau, taussa
  # Scalars to work around OpenACC/OMP issues
  # real(wp) :: minSize,  maxSize

  # ----------------------------------------
  #
  # Error checking
  #
  # ----------------------------------------
  allocated(this.aero_dust_tbl) || return "aerosol optics: no data has been initialized"

  ncol = size(aero_type,1)
  nlay = size(aero_type,2)
  npair= size(this.merra_aero_bin_lims,1)
  nbin = size(this.merra_aero_bin_lims,2)
  nrh  = size(this.aero_rh,1)
  nval = size(this.aero_dust_tbl,1)
  nbnd = size(this.aero_dust_tbl,3)

  # !$acc        update host(this.merra_aero_bin_lims)
  # !$omp target update from(this.merra_aero_bin_lims)
  minSize = this.merra_aero_bin_lims[1,1]
  maxSize = this.merra_aero_bin_lims[2,nbin]

  #
  # Array sizes
  #
  if check_extents
    extents_are(aero_type, ncol, nlay) || return "aerosol optics: aero_type isn't consistenly sized"
    extents_are(aero_size, ncol, nlay) || return "aerosol optics: aero_size isn't consistenly sized"
    extents_are(aero_mass, ncol, nlay) || return "aerosol optics: aero_mass isn't consistenly sized"
    extents_are(relhum, ncol, nlay) || return "aerosol optics: relhum isn't consistenly sized"
    if get_ncol(optical_props) ≠ ncol || get_nlay(optical_props) ≠ nlay
      return "aerosol optics: optical_props have wrong extents"
    end
  end

  #
  # Spectral consistency
  #
  if check_values
    bands_are_equal(this, optical_props) || return "aerosol optics: optical properties don't have the same band structure"
    get_nband(optical_props) ≠ get_ngpt(optical_props) && return "aerosol optics: optical properties must be requested by band not g-points"
    any_int_vals_outside_2D(aero_type, merra_aero_none, merra_ntype) && return "aerosol optics: aerosol type is out of bounds"
  end

  # !$acc data        copyin(aero_type, aero_size, aero_mass, relhum)
  # !$omp target data map(to:aero_type, aero_size, aero_mass, relhum)
  #
  # Aerosol mask; don't need aerosol optics if there's no aerosol
  #
  aeromsk = zeros(size(aero_type))
  # !$acc data           create(aeromsk)
  # !$omp target data map(alloc:aeromsk)
  # !$acc              parallel loop default(present) collapse(2)
  # !$omp target teams distribute parallel do simd collapse(2)
  for ilay in 1:nlay
    for icol in 1:ncol
      aeromsk[icol,ilay] = aero_type[icol,ilay] > 0
    end
  end

  #
  # Aerosol size, relative humidity
  #
  if check_values
    any_vals_outside(aero_size, aeromsk, minSize, maxSize) && return "aerosol optics: requested aerosol size is out of bounds"
    any_vals_outside(relhum,    aeromsk, 0, 1) && return "aerosol optics: relative humidity fraction is out of bounds"
  end
  # Release aerosol mask
  # !$acc end data
  # !$omp end target data

  error_msg ≠ "" && return error_msg

  # real(wp),    dimension(size(aero_type,1), size(aero_type,2), size(this.aero_dust_tbl,3)) :: &
  #              atau, ataussa, ataussag

    atau = zeros((size(aero_type)..., size(this.aero_dust_tbl,3)))
    ataussa = zeros((size(aero_type)..., size(this.aero_dust_tbl,3)))
    ataussag = zeros((size(aero_type)..., size(this.aero_dust_tbl,3)))
    # !$acc data           create(atau, ataussa, ataussag)
    # !$omp target data map(alloc:atau, ataussa, ataussag)
    # ----------------------------------------
    #
    # The lookup tables determining extinction coefficient, single-scattering albedo,
    #   and asymmetry parameter g as a function of aerosol type, aerosol size and
    #   relative humidity.
    # We compute the optical depth tau (= exintinction coeff * aerosol mass ) and the
    #    products tau*ssa and tau*ssa*g separately for each aerosol type requested.
    # These are used to determine the aerosol optical properties.
    if allocated(this.aero_dust_tbl)
      #
      # Aerosol
      #
       compute_all_from_table(ncol, nlay, npair, nval, nrh, nbin, nbnd,
                              aero_type, aero_size, aero_mass, relhum,
                              this.merra_aero_bin_lims, this.aero_rh,
                              this.aero_dust_tbl, this.aero_salt_tbl, this.aero_sulf_tbl,
                              this.aero_bcar_rh_tbl, this.aero_bcar_tbl,
                              this.aero_ocar_rh_tbl, this.aero_ocar_tbl,
                              atau, ataussa, ataussag)

    end

    #
    # Derive total aerosol optical properties
    #   See also the increment routines in mo_optical_props_kernels
    #
    if optical_props isa ty_optical_props_1scl
      # !$acc parallel loop gang vector default(present) collapse(3) &
      # !$acc               copyin(optical_props) copyout(optical_props.tau)
      # !$omp target teams distribute parallel do simd collapse(3) &
      # !$omp map(from:optical_props.tau)
      for ibnd in 1:nbnd
        for ilay in 1:nlay
          for icol in 1:ncol
            # Absorption optical depth  = (1-ssa) * tau = tau - taussa
            optical_props.tau[icol,ilay,ibnd] = (atau[icol,ilay,ibnd] - ataussa[icol,ilay,ibnd])
          end
        end
      end
    elseif optical_props isa ty_optical_props_2str
      # !$acc parallel loop gang vector default(present) collapse(3) &
      # !$acc               copyin(optical_props) copyout(optical_props.tau, optical_props.ssa, optical_props.g)
      # !$omp target teams distribute parallel do simd collapse(3) &
      # !$omp map(from:optical_props.tau, optical_props.ssa, optical_props.g)
      for ibnd in 1:nbnd
        for ilay in 1:nlay
          for icol in 1:ncol
            tau    = atau[icol,ilay,ibnd]
            taussa = ataussa[icol,ilay,ibnd]
            optical_props.tau[icol,ilay,ibnd] = tau
            optical_props.ssa[icol,ilay,ibnd] = taussa / max(fepsilon(tau), tau)
            optical_props.g[icol,ilay,ibnd] = (ataussag[icol,ilay,ibnd]) / max(fepsilon(tau), taussa)
          end
        end
      end
    elseif optical_props isa ty_optical_props_nstr
    else
      return "aerosol optics: n-stream calculations not yet supported"
    end
    # !$acc end data
    # !$omp end target data
  # !$acc end data
  # !$omp end target data
end
#--------------------------------------------------------------------------------------------------------------------
#
# Ancillary functions
#
#--------------------------------------------------------------------------------------------------------------------
#
# For asize dimension, select asize bin appropriate for the requested aerosol size.
# For rh dimension, linearly interpolate values from a lookup table with "nrh"
#   unevenly-spaced elements "aero_rh". The last dimension for all tables is band.
# Returns zero where no aerosol is present.
#
function compute_all_from_table(ncol, nlay, npair, nval, nrh, nbin, nbnd,
    type, asize, mass, rh,
    merra_aero_bin_lims, aero_rh,
    aero_dust_tbl, aero_salt_tbl, aero_sulf_tbl,
    aero_bcar_rh_tbl, aero_bcar_tbl,
    aero_ocar_rh_tbl, aero_ocar_tbl,
    tau, taussa, taussag)

  FT = eltype(tau)
  # integer,                               intent(in) :: ncol, nlay, npair, nval, nrh, nbin, nbnd
  # integer,     dimension(ncol,nlay),     intent(in) :: type
  # real(wp),    dimension(ncol,nlay),     intent(in) :: asize, mass, rh

  # real(wp),    dimension(npair,nbin),    intent(in) :: merra_aero_bin_lims
  # real(wp),    dimension(nrh),           intent(in) :: aero_rh

  # real(wp),    dimension(nval,    nbin,nbnd), intent(in) :: aero_dust_tbl
  # real(wp),    dimension(nrh,nval,nbin,nbnd), intent(in) :: aero_salt_tbl
  # real(wp),    dimension(nrh,nval,     nbnd), intent(in) :: aero_sulf_tbl
  # real(wp),    dimension(nrh,nval,     nbnd), intent(in) :: aero_bcar_rh_tbl
  # real(wp),    dimension(nval,         nbnd), intent(in) :: aero_bcar_tbl
  # real(wp),    dimension(nrh,nval,     nbnd), intent(in) :: aero_ocar_rh_tbl
  # real(wp),    dimension(nval,         nbnd), intent(in) :: aero_ocar_tbl

  # real(wp),    dimension(ncol,nlay,nbnd), intent(out) :: tau, taussa, taussag
  # ! ---------------------------
  # integer  :: icol, ilay, ibnd, ibin, i
  # integer  :: itype, irh1, irh2
  # real(wp) :: drh0, drh1, rdrh
  # real(wp) :: t, ts, tsg  ! tau, tau*ssa, tau*ssa*g
  # ! ---------------------------
  # !$acc parallel loop gang vector default(present) collapse(3)
  # !$omp target teams distribute parallel do simd collapse(3)
  ibin = -1
  for ibnd in 1:nbnd
    for ilay in 1:nlay
      for icol in 1:ncol
        # Sequential loop to find asize bin
        for i in 1:nbin
           if (asize[icol,ilay] ≥ merra_aero_bin_lims[1,i] &&
               asize[icol,ilay] ≤ merra_aero_bin_lims[2,i])
              ibin = i
           end
        end
        itype = type[icol,ilay]
        # relative humidity linear interpolation coefficients
        if itype ≠ merra_aero_none
           irh2 = 1
           while rh[icol,ilay] > aero_rh[irh2]
              irh2 = irh2 + 1
              (irh2 > nrh) && break
           end
           irh1 = max(1, irh2-1)
           irh2 = min(nrh, irh2)
           drh0 = aero_rh[irh2] - aero_rh[irh1]
           drh1 = rh[icol,ilay] - aero_rh[irh1]
           if (irh1 == irh2)
              rdrh = FT(0)
           else
              rdrh = drh1 / drh0
           end
        end

        # Set aerosol optical properties where aerosol present. Use aerosol type array as the mask.
         # dust
         if itype == merra_aero_dust
           tau[icol,ilay,ibnd]     = mass[icol,ilay]        * aero_dust_tbl[ext,ibin,ibnd]
           taussa[icol,ilay,ibnd]  = tau[icol,ilay,ibnd]    * aero_dust_tbl[ssa,ibin,ibnd]
           taussag[icol,ilay,ibnd] = taussa[icol,ilay,ibnd] * aero_dust_tbl[g,ibin,ibnd]
         # sea-salt
         elseif itype == merra_aero_salt
           tau[icol,ilay,ibnd]     = mass[icol,ilay]        * linear_interp_aero_table(aero_salt_tbl[:,ext,ibin,ibnd],irh1,irh2,rdrh)
           taussa[icol,ilay,ibnd]  = tau[icol,ilay,ibnd]    * linear_interp_aero_table(aero_salt_tbl[:,ssa,ibin,ibnd],irh1,irh2,rdrh)
           taussag[icol,ilay,ibnd] = taussa[icol,ilay,ibnd] * linear_interp_aero_table(aero_salt_tbl[:,g,  ibin,ibnd],irh1,irh2,rdrh)

         # sulfate
         elseif itype == merra_aero_sulf
           tau[icol,ilay,ibnd]     = mass[icol,ilay]        * linear_interp_aero_table(aero_sulf_tbl[:,ext,ibnd],irh1,irh2,rdrh)
           taussa[icol,ilay,ibnd]  = tau[icol,ilay,ibnd]    * linear_interp_aero_table(aero_sulf_tbl[:,ssa,ibnd],irh1,irh2,rdrh)
           taussag[icol,ilay,ibnd] = taussa[icol,ilay,ibnd] * linear_interp_aero_table(aero_sulf_tbl[:,g,  ibnd],irh1,irh2,rdrh)
         # black carbon - hydrophilic
         elseif itype == merra_aero_bcar_rh
           tau[icol,ilay,ibnd]     = mass[icol,ilay]        * linear_interp_aero_table(aero_bcar_rh_tbl[:,ext,ibnd],irh1,irh2,rdrh)
           taussa[icol,ilay,ibnd]  = tau[icol,ilay,ibnd]    * linear_interp_aero_table(aero_bcar_rh_tbl[:,ssa,ibnd],irh1,irh2,rdrh)
           taussag[icol,ilay,ibnd] = taussa[icol,ilay,ibnd] * linear_interp_aero_table(aero_bcar_rh_tbl[:,g,  ibnd],irh1,irh2,rdrh)
         # black carbon - hydrophobic
         elseif itype == merra_aero_bcar
           tau[icol,ilay,ibnd]     = mass[icol,ilay]        * aero_bcar_tbl[ext,ibnd]
           taussa[icol,ilay,ibnd]  = tau[icol,ilay,ibnd]    * aero_bcar_tbl[ssa,ibnd]
           taussag[icol,ilay,ibnd] = taussa[icol,ilay,ibnd] * aero_bcar_tbl[g,  ibnd]
         # organic carbon - hydrophilic
         elseif itype == merra_aero_ocar_rh
           tau[icol,ilay,ibnd]     = mass[icol,ilay]        * linear_interp_aero_table(aero_ocar_rh_tbl[:,ext,ibnd],irh1,irh2,rdrh)
           taussa[icol,ilay,ibnd]  = tau[icol,ilay,ibnd]    * linear_interp_aero_table(aero_ocar_rh_tbl[:,ssa,ibnd],irh1,irh2,rdrh)
           taussag[icol,ilay,ibnd] = taussa[icol,ilay,ibnd] * linear_interp_aero_table(aero_ocar_rh_tbl[:,g,  ibnd],irh1,irh2,rdrh)
         # organic carbon - hydrophobic
         elseif itype == merra_aero_ocar
           tau[icol,ilay,ibnd]     = mass[icol,ilay]        * aero_ocar_tbl[ext,ibnd]
           taussa[icol,ilay,ibnd]  = tau[icol,ilay,ibnd]    * aero_ocar_tbl[ssa,ibnd]
           taussag[icol,ilay,ibnd] = taussa[icol,ilay,ibnd] * aero_ocar_tbl[g,  ibnd]
         # no aerosol
         else
           tau[icol,ilay,ibnd]     = 0
           taussa[icol,ilay,ibnd]  = 0
           taussag[icol,ilay,ibnd] = 0

        end

      end
    end
  end
end
#--------------------------------------------------------------------------------------------------------------------
#
# Function for linearly interpolating MERRA aerosol optics tables in the rh dimension for
# a single parameter, aerosol type, spectral band, and size bin. Interpolation is performed
# only where aerosol in present using aerosol type as the mask.
#
function linear_interp_aero_table(table::AbstractVector, index1::Int, index2::Int, weight::Real)
  # !$acc routine seq
  # !$omp declare target

  # integer,                intent(in) :: index1, index2
  # real(wp),               intent(in) :: weight
  # real(wp), dimension(:), intent(in) :: table

  # real(wp) :: value

  return table[index1] + weight * (table(index2) - table(index1))

end
# ----------------------------------------------------------
function any_int_vals_outside_2D(array::AbstractArray, checkMin::Int, checkMax::Int)::Bool
  # integer, dimension(:,:), intent(in) :: array
  # integer,                 intent(in) :: checkMin, checkMax

  # integer :: minValue, maxValue

  # !$acc kernels copyin(array)
  # !$omp target map(to:array) map(from:minValue, maxValue)
  minValue = minval(array)
  maxValue = maxval(array)
  # !$acc end kernels
  # !$omp end target
  any_int_vals_outside_2D = minValue < checkMin || maxValue > checkMax

end
