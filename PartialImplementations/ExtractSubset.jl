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
    for igpt in 1:ngpt
        for ilay in 1:nlay
            for icol in colS:colE
                array_out[icol - colS + 1, ilay, igpt] = array_in[icol, ilay, igpt]
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

    for igpt in 1:ngpt
        for ilay in 1:nlay
            for icol in colS:colE
                for imom in 1:nmom
                    array_out[imom, icol - colS + 1, ilay, igpt] = array_in[imom, icol, ilay, igpt]
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

    for igpt in 1:ngpt
        for ilay in 1:nlay
            for icol in colS:colE
                tau_out[icol - colS + 1, ilay, igpt] = tau_in[icol, ilay, igpt] * (FT(1) - ssa_in[icol, ilay, igpt])
            end
        end
    end
end


# Extract a subset of n columns starting with column 'start'
function get_subset_range(this::GasConcs, start::I, n::I, subset::GasConcs) where {I}
    # class(GasConcs),      intent(in   ) :: this
    # integer,                  intent(in   ) :: start, n
    # class(GasConcs),      intent(inout) :: subset
    @assert !(n <= 0)
    @assert !(start < 1)
    @assert !(this.ncol ≠ nothing && start > this.ncol || start + n - 1 > this.ncol)

    this.nlay = 0
    this.ncol = 0
    # reset!(subset)
    # These two arrays should be the same length
    subset.gas_name = Array(undef, size(this.gas_name))
    subset.concs = Array(undef, size(this.concs))
    subset.nlay = this.nlay
    subset.ncol = this.ncol ≠ nothing ? n : 0
    subset.gas_name[:] .= this.gas_name[:]

    for i in 1:size(this.gas_name)
        #
        # Preserve scalar/1D/2D representation in subset,
        #   but need to ensure at least extent 1 in col dimension (ncol = 0 means no gas exploits this dimension)
        #
        allocate(
            subset.concs[i].conc(
                min(max(subset.ncol, 1), size(this.concs[i].conc, 1)),
                min(subset.nlay, size(this.concs[i].conc, 2)),
            ),
        )
        if size(this.concs[i].conc, 1) > 1      # Concentration stored as 2D
            subset.concs[i].conc[:, :] .= this.concs[i].conc[start:(start + n - 1), :]
        else
            subset.concs[i].conc[:, :] .= this.concs[i].conc[:, :]
        end
    end

end
