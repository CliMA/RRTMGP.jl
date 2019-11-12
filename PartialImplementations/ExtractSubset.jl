#### Extract subset
# -------------------------------------------------------------------------------------------------
#
# Subsetting, meaning extracting some portion of the 3D domain
#
# -------------------------------------------------------------------------------------------------
function extract_subset_dim1_3d(ncol, nlay, ngpt, array_in, colS, colE, array_out)#
#    bind (C, name="extract_subset_dim1_3d")
#    integer,                             intent(in ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol,nlay,ngpt), intent(in ) :: array_in
#    integer,                             intent(in ) :: colS, colE
#    real(FT), dimension(colE-colS+1,
#                             nlay,ngpt), intent(out) :: array_out

#    integer :: icol, ilay, igpt
  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = colS:colE
        array_out[icol-colS+1, ilay, igpt] = array_in[icol, ilay, igpt]
      end
    end
  end

end
# ---------------------------------
function extract_subset_dim2_4d(nmom, ncol, nlay, ngpt, array_in, colS, colE, array_out)#
#    bind (C, name="extract_subset_dim2_4d")
#    integer,                                  intent(in ) :: nmom, ncol, nlay, ngpt
#    real(FT), dimension(nmom,ncol,nlay,ngpt), intent(in ) :: array_in
#    integer,                                  intent(in ) :: colS, colE
#    real(FT), dimension(nmom,colE-colS+1,
#                                  nlay,ngpt), intent(out) :: array_out

#    integer :: icol, ilay, igpt, imom

  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = colS:colE
        for imom = 1:nmom
          array_out[imom, icol-colS+1, ilay, igpt] = array_in[imom, icol, ilay, igpt]
        end
      end
    end
  end

end

#
# Extract the absorption optical thickness which requires mulitplying by 1 - ssa
#
#    integer,                             intent(in ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol,nlay,ngpt), intent(in ) :: tau_in, ssa_in
#    integer,                             intent(in ) :: colS, colE
#    real(FT), dimension(colE-colS+1,
#                             nlay,ngpt), intent(out) :: tau_out
#    integer :: icol, ilay, igpt
function extract_subset_absorption_tau(ncol, nlay, ngpt, tau_in, ssa_in, colS, colE, tau_out)
  FT = eltype(tau_in)

  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = colS:colE
        tau_out[icol-colS+1, ilay, igpt] = tau_in[icol, ilay, igpt] * ( FT(1) - ssa_in[icol, ilay, igpt] )
      end
    end
  end
end


# Extract a subset of n columns starting with column 'start'
function get_subset_range(this::ty_gas_concs, start::I, n::I, subset::ty_gas_concs) where I
  # class(ty_gas_concs),      intent(in   ) :: this
  # integer,                  intent(in   ) :: start, n
  # class(ty_gas_concs),      intent(inout) :: subset
  @assert !(n <= 0)
  @assert !(start < 1 )
  @assert !(this.ncol ≠ nothing && start > this.ncol || start+n-1 > this.ncol )

  this.nlay = 0
  this.ncol = 0
  # reset!(subset)
  # These two arrays should be the same length
  subset.gas_name = Array(undef, size(this.gas_name))
  subset.concs = Array(undef, size(this.concs))
  subset.nlay = this.nlay
  subset.ncol = fmerge(n, 0, this.ncol ≠ nothing)
  subset.gas_name[:] .= this.gas_name[:]

  for i = 1:size(this.gas_name)
    #
    # Preserve scalar/1D/2D representation in subset,
    #   but need to ensure at least extent 1 in col dimension (ncol = 0 means no gas exploits this dimension)
    #
    allocate(subset.concs[i].conc(min(max(subset.ncol,1), size(this.concs[i].conc, 1)),
                                        min(    subset.nlay,    size(this.concs[i].conc, 2))))
    if size(this.concs[i].conc, 1) > 1      # Concentration stored as 2D
      subset.concs[i].conc[:,:] .= this.concs[i].conc[start:(start+n-1),:]
    else
      subset.concs[i].conc[:,:] .= this.concs[i].conc[:,:]
    end
  end

end
