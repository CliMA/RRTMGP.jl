# This file provides a simple implementation of sampling for the
# Monte Carlo Independent Pixel Approximation (McICA, doi:10.1029/2002jd003322)
# Cloud optical properties, defined by band and assumed homogenous within each cell (column/layer),
# are randomly sampled to preserve the mean cloud fraction and one of several possible overlap assumptions
# Users supply random numbers with order ngpt,nlay,ncol
# These are only accessed if cloud_fraction(icol,ilay) > 0 so many values don't need to be filled in

"""
    draw_samples

Apply a T/F sampled cloud mask to cloud optical properties defined by band to produce
   McICA-sampled cloud optical properties
"""
function draw_samples(cloud_mask,clouds,clouds_sampled)
  # logical, dimension(:,:,:),      intent(in   ) :: cloud_mask     # Dimensions ncol,nlay,ngpt
  # class(ty_optical_props_arry),   intent(in   ) :: clouds         # Defined by band
  # class(ty_optical_props_arry),   intent(inout) :: clouds_sampled # Defined by g-point
  # character(len=128)                            :: error_msg
  # # ------------------------
  # integer :: ncol,nlay,nbnd,ngpt
  # integer :: imom
  # # ------------------------

  #
  # Variables clouds and clouds_sampled have to be of the same type (have the same set of fields)
  #   nstr isn't supported
  #   2str is checked at assignment
  #

  if clouds isa ty_optical_props_1scl
    if clouds_sampled isa ty_optical_props_2str
      error("draw_samples: by-band and sampled cloud properties need to be the same variable type")
    end
  end

  #
  # Spectral discretization
  #
  @assert bands_are_equal(clouds, clouds_sampled)

  #
  # Array extents
  #
  ncol = get_ncol(clouds)
  nlay = get_nlay(clouds)
  nbnd = get_nband(clouds)
  ngpt = get_ngpt(clouds_sampled)
  if any([size(cloud_mask,1), size(cloud_mask,2), size(cloud_mask,3)] ≠ [ncol,nlay,ngpt])
    error("draw_samples: cloud mask and cloud optical properties have different ncol and/or nlay")
  end
  if (any([get_ncol(clouds_sampled), get_nlay(clouds_sampled)] ≠ [ncol,nlay]))
    error("draw_samples: sampled/unsampled cloud optical properties have different ncol and/or nlay")
  end
  # ------------------------
  #
  # Finally - sample fields according to the cloud mask
  #
  # Optical depth assignment works for 1scl, 2str (also nstr)
  apply_cloud_mask!(ncol,nlay,nbnd,ngpt,get_band_lims_gpoint(clouds_sampled),cloud_mask,clouds.tau,clouds_sampled.tau)
  #
  # For 2-stream
  #
  if clouds isa ty_optical_props_2str
    if clouds_sampled isa ty_optical_props_2str
      apply_cloud_mask!(ncol,nlay,nbnd,ngpt,get_band_lims_gpoint(clouds_sampled),cloud_mask,clouds.ssa,clouds_sampled.ssa)
      apply_cloud_mask!(ncol,nlay,nbnd,ngpt,get_band_lims_gpoint(clouds_sampled),cloud_mask,clouds.g,  clouds_sampled.g  )
    else
      error("draw_samples: by-band and sampled cloud properties need to be the same variable type")
    end
  end
end


"""
    sampled_mask_max_ran!

Generate a McICA-sampled cloud mask for maximum-random overlap
"""
function sampled_mask_max_ran!(randoms,cloud_frac,cloud_mask)
#   real(FT), dimension(:,:,:),    intent(in ) :: randoms    #ngpt,nlay,ncol
#   real(FT), dimension(:,:),      intent(in ) :: cloud_frac # ncol,nlay
#   logical,  dimension(:,:,:),    intent(out) :: cloud_mask # ncol,nlay,ngpt
#   character(len=128)                         :: error_msg
# # ------------------------
#   integer                              :: ncol, nlay, ngpt, icol, ilay, igpt
#   integer                              :: cloud_lay_fst, cloud_lay_lst
#   real(FT), dimension(size(randoms,1)) :: local_rands
#   logical,  dimension(size(randoms,2)) :: cloud_mask_layer

  # ------------------------
  #
  # Error checking
  #
  FT = eltype(randoms)
  ncol = size(randoms, 3)
  nlay = size(randoms, 2)
  ngpt = size(randoms, 1)
  if any([ncol,nlay] ≠ [size(cloud_frac, 1),size(cloud_frac, 2)])
    error("sampled_mask_max_ran: sizes of randoms(ngpt,nlay,ncol) and cloud_frac(ncol,nlay) are inconsistent")
  end
  if any([ncol,nlay,ngpt] ≠ [size(cloud_mask, 1),size(cloud_mask, 2), size(cloud_mask,3)])
    error("sampled_mask_max_ran: sizes of randoms(ngpt,nlay,ncol) and cloud_mask(ncol,nlay,ngpt) are inconsistent")
  end
  if any(cloud_frac > FT(1)) || any(cloud_frac < FT(0))
    error("sampled_mask_max_ran: cloud fraction values out of range [0,1]")
  end
  #
  # We chould check the random numbers but that would be computationally heavy
  #
  # ------------------------
  #
  # Construct the cloud mask for each column
  #
  for icol = 1:ncol
    cloud_mask_layer[1:nlay] = cloud_frac[icol,1:nlay] > FT(0)
    if(!any(cloud_mask_layer))
      cloud_mask[icol,1:nlay,1:ngpt] = false
      continue
    end
    cloud_lay_fst = findloc(cloud_mask_layer, true, dim=1)
    cloud_lay_lst = findloc(cloud_mask_layer, true, dim=1, back = true)
    cloud_mask[icol,1:cloud_lay_fst,1:ngpt] = false

    ilay = cloud_lay_fst
    local_rands[1:ngpt] = randoms[1:ngpt,cloud_lay_fst,icol]
    cloud_mask[icol,ilay,1:ngpt] = local_rands[1:ngpt] > (FT(1) - cloud_frac[icol,ilay])
    for ilay = cloud_lay_fst+1:cloud_lay_lst
      if cloud_mask_layer[ilay]
        #
        # Max-random overlap:
        #   new  random deviates if the adjacent layer isn't cloudy
        #   same random deviates if the adjacent layer is    cloudy
        #
        if(!cloud_mask_layer[ilay-1])
          local_rands[1:ngpt] = randoms[1:ngpt,ilay,icol]
        end
        cloud_mask[icol,ilay,1:ngpt] = local_rands[1:ngpt] > (FT(1) - cloud_frac[icol,ilay])
      else
        cloud_mask[icol,ilay,1:ngpt] = false
      end
    end

    cloud_mask[icol,cloud_lay_lst+1:nlay,1:ngpt] = false
  end

end

"""
    sampled_mask_exp_ran

 Generate a McICA-sampled cloud mask for exponential-random overlap
   The overlap parameter alpha is defined between pairs of layers
   for layer i, alpha(i) describes the overlap betwen cloud_frac(i) and cloud_frac(i+1)
   By skipping layers with 0 cloud fraction the code forces alpha(i) = 0 for cloud_frac(i) = 0.
"""
function sampled_mask_exp_ran(randoms,cloud_frac,overlap_param,cloud_mask)
  # real(FT), dimension(:,:,:), intent(in ) :: randoms       # ngpt,nlay,ncol
  # real(FT), dimension(:,:),   intent(in ) :: cloud_frac    # ncol,nlay
  # real(FT), dimension(:,:),   intent(in ) :: overlap_param # ncol,nlay-1
  # logical,  dimension(:,:,:), intent(out) :: cloud_mask    # ncol,nlay,ngpt
  # character(len=128)                      :: error_msg
  # # ------------------------
  # integer                              :: ncol, nlay, ngpt, icol, ilay, igpt
  # integer                              :: cloud_lay_fst, cloud_lay_lst
  # real(FT)                             :: rho # correlation coefficient
  # real(FT), dimension(size(randoms,1)) :: local_rands
  # logical,  dimension(size(randoms,2)) :: cloud_mask_layer
  # ------------------------
  #
  # Error checking
  #
  ncol = size(randoms, 3)
  nlay = size(randoms, 2)
  ngpt = size(randoms, 1)
  if(any([ncol,nlay] ≠ [size(cloud_frac, 1),size(cloud_frac, 2)]))
    error("sampled_mask_max_ran: sizes of randoms(ngpt,nlay,ncol) and cloud_frac(ncol,nlay) are inconsistent")
  end
  if(any([ncol,nlay-1] ≠ [size(overlap_param, 1),size(overlap_param, 2)]))
    error("sampled_mask_max_ran: sizes of randoms(ngpt,nlay,ncol) and overlap_param(ncol,nlay-1) are inconsistent")
  end
  if(any([ncol,nlay,ngpt] ≠ [size(cloud_mask, 1),size(cloud_mask, 2), size(cloud_mask,3)]))
    error("sampled_mask_max_ran: sizes of randoms(ngpt,nlay,ncol) and cloud_mask(ncol,nlay,ngpt) are inconsistent")
  end

  if(any(cloud_frac > FT(1)) || any(cloud_frac < FT(0)))
    error("sampled_mask_max_ran: cloud fraction values out of range [0,1]")
  end
  if(any(overlap_param > FT(1)) || any(overlap_param < -FT(1)))
    error("sampled_mask_max_ran: overlap_param values out of range [-1,1]")
  end
  #
  # We chould check the random numbers but that would be computationally heavy
  #
  # ------------------------
  # Construct the cloud mask for each column
  #
  for icol = 1:ncol
    cloud_mask_layer[1:nlay] .= cloud_frac[icol,1:nlay] .> FT(0)
    if !any(cloud_mask_layer)
      cloud_mask[icol,1:nlay,1:ngpt] .= false
      continue
    end
    cloud_lay_fst = findloc(cloud_mask_layer, true, dim=1)
    cloud_lay_lst = findloc(cloud_mask_layer, true, dim=1, back = true)
    cloud_mask[icol,1:cloud_lay_fst,1:ngpt] = false

    ilay = cloud_lay_fst
    local_rands[1:ngpt] = randoms[1:ngpt,ilay,icol]
    cloud_mask[icol,ilay,1:ngpt] = local_rands[1:ngpt] > (FT(1) - cloud_frac[icol,ilay])
    for ilay = cloud_lay_fst+1:cloud_lay_lst
      if cloud_mask_layer[ilay]
        #
        # Exponential-random overlap:
        #   new  random deviates if the adjacent layer isn't cloudy
        #   correlated  deviates if the adjacent layer is    cloudy
        #
        if cloud_mask_layer[ilay-1]
          #
          # Create random deviates correlated between this layer and the previous layer
          #    (have to remove mean value before enforcing correlation)
          #
          rho = overlap_param[icol,ilay-1]
          local_rands[1:ngpt] =  rho*(local_rands[1:ngpt]      -FT(0.5)) +
                 sqrt(FT(1)-rho*rho)*(randoms[1:ngpt,ilay,icol]-FT(0.5)) + FT(0.5)
        else
          local_rands[1:ngpt] = randoms[1:ngpt,ilay,icol]
        end
        cloud_mask[icol,ilay,1:ngpt] = local_rands[1:ngpt] > (FT(1) - cloud_frac[icol,ilay])
      end
    end

    cloud_mask[icol,cloud_lay_lst+1:nlay, 1:ngpt] = false
  end

end

"""
    apply_cloud_mask

Apply a true/false cloud mask to a homogeneous field
"""
function apply_cloud_mask(ncol,nlay,nbnd,ngpt,band_lims_gpt,cloud_mask,input_field,sampled_field)
  # integer,                                intent(in ) :: ncol,nlay,nbnd,ngpt
  # integer,     dimension(2,nbnd),         intent(in ) :: band_lims_gpt
  # logical,     dimension(ncol,nlay,ngpt), intent(in ) :: cloud_mask
  # real(FT),    dimension(ncol,nlay,nbnd), intent(in ) :: input_field
  # real(FT),    dimension(ncol,nlay,ngpt), intent(out) :: sampled_field

  # integer :: icol,ilay,ibnd,igpt

  for ibnd = 1:nbnd
    for igpt = band_lims_gpt[1,ibnd]:band_lims_gpt[2,ibnd]
      for ilay = 1:nlay
        sampled_field[1:ncol,ilay,igpt] = merge(input_field[1:ncol,ilay,ibnd], FT(0), cloud_mask[1:ncol,ilay,igpt])
      end
    end
  end
end
