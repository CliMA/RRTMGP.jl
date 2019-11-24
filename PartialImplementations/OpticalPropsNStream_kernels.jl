#### n-stream kernels
#  - delta-scaling
#  - adding two sets of properties
#  - extracting subsets
#  - validity checking

# ---------------------------------
# increment 1scalar by nstream
function increment_by_gpoint!(op_1::OneScalar{FT}, op_2::NStream{FT}) where FT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau
#    real(FT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa

#    integer  :: icol, ilay, igpt

  op_1.tau = op_1.tau .+ op_2.tau * ( FT(1) - op_2.ssa)
end

# ---------------------------------
# increment 2stream by nstream
function increment_by_gpoint!(op_1::TwoStream{FT}, op_2::NStream{FT}) where FT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom2
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa, op_1.g
#    real(FT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa
#    real(FT), dimension(nmom2,
#                        ncol,nlay,ngpt), intent(in   ) :: op_2.p

#    integer  :: icol, ilay, igpt
#    real(FT) :: tau12, tauscat12

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)

  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = 1:ncol
        # t=op_1.tau + op_2.tau
        tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,igpt]
        # w=(op_1.tau*op_1.ssa + op_2.tau*op_2.ssa) / t
        tauscat12 =
           op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
           op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt]
        op_1.g[icol,ilay,igpt] =
          (op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.g[   icol,ilay,igpt]+
           op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt] * op_2.p[1, icol,ilay,igpt]) / max(ε_machine(FT),tauscat12)
        op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(FT),tau12)
        op_1.tau[icol,ilay,igpt] = tau12
      end
    end
  end
end

# increment nstream by 1scalar
function increment_by_gpoint!(op_1::NStream{FT}, op_2::OneScalar{FT}) where FT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(FT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau

#    integer  :: icol, ilay, igpt
#    real(FT) :: tau12

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)

  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = 1:ncol
        tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,igpt]
        op_1.ssa[icol,ilay,igpt] = op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(ε_machine(FT),tau12)
        op_1.tau[icol,ilay,igpt] = tau12
        # p is unchanged
      end
    end
  end
end

# ---------------------------------
# increment nstream by 2stream
function increment_by_gpoint!(op_1::NStream{FT}, op_2::TwoStream{FT}) where FT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom1
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(FT), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: op_1.p
#    real(FT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa, op_2.g

#    integer  :: icol, ilay, igpt
#    real(FT) :: tau12, tauscat12
#    real(FT), dimension(nmom1) :: temp_moms # TK
#    integer  :: imom  #TK

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)

  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = 1:ncol
        tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,igpt]
        tauscat12 =
           op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
           op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt]
        #
        # Here assume Henyey-Greenstein
        #
        temp_moms[1] = op_2.g[icol,ilay,igpt]
        for imom = 2:nmom1
          temp_moms[imom] = temp_moms[imom-1] * op_2.g[icol,ilay,igpt]
        end
        op_1.p[1:nmom1, icol,ilay,igpt] =
            (op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.p[1:nmom1, icol,ilay,igpt] +
             op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt] * temp_moms[1:nmom1]  ) / max(ε_machine(FT),tauscat12)
        op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(FT),tau12)
        op_1.tau[icol,ilay,igpt] = tau12
      end
    end
  end
end


# ---------------------------------
# increment nstream by nstream
function increment_by_gpoint!(op_1::NStream{FT}, op_2::NStream{FT}) where FT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom1, nmom2
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(FT), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: op_1.p
#    real(FT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa
#    real(FT), dimension(nmom2,
#                        ncol,nlay,ngpt), intent(in   ) :: op_2.p

#    integer  :: icol, ilay, igpt, mom_lim
#    real(FT) :: tau12, tauscat12

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)

  mom_lim = min(nmom1, nmom2)
  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = 1:ncol
        tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,igpt]
        tauscat12 =
           op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
           op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt]
        #
        # If op2 has more moments than op1 these are ignored;
        #   if it has fewer moments the higher orders are assumed to be 0
        #
        op_1.p[1:mom_lim, icol,ilay,igpt] =
            (op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.p[1:mom_lim, icol,ilay,igpt] +
             op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt] * op_2.p[1:mom_lim, icol,ilay,igpt]) / max( ε_machine(FT),tauscat12)
        op_1.ssa[icol,ilay,igpt] = tauscat12 / max(eps,tau12)
        op_1.tau[icol,ilay,igpt] = tau12
      end
    end
  end
end

# increment 1scalar by nstream
function increment_bybnd!(op_1::OneScalar{FT}, op_2::NStream{FT}) where FT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)
  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      op_1.tau[:,:,igpt] .= op_1.tau[:,:,igpt] .+ op_2.tau[:,:,ibnd] .* ( FT(1) .- op_2.ssa[:,:,ibnd])
    end
  end
end

"""
    increment_bybnd!(op_1::TwoStream{FT}, op_2::NStream{FT}) where FT

integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom2, nbnd
real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa, op_1.g
real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: op_2.tau, op_2.ssa
real(FT), dimension(nmom2,
                    ncol,nlay,nbnd), intent(in   ) :: op_2.p
integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

integer  :: icol, ilay, igpt, ibnd
real(FT) :: tau12, tauscat12
"""
function increment_bybnd!(op_1::TwoStream{FT}, op_2::NStream{FT}) where FT

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)
  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      for ilay = 1:nlay
        for icol = 1:ncol
          # t=op_1.tau + op_2.tau
          tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,ibnd]
          # w=(op_1.tau*op_1.ssa + op_2.tau*op_2.ssa) / t
          tauscat12 =
             op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
             op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd]
          op_1.g[icol,ilay,igpt] =
            (op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.g[   icol,ilay,igpt]+
             op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd] * op_2.p[1, icol,ilay,ibnd]) / max(ε_machine(FT),tauscat12)
          op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(FT),tau12)
          op_1.tau[icol,ilay,igpt] = tau12
        end
      end
    end
  end
end

# increment nstream by 1scalar
function increment_bybnd!(op_1::NStream{FT}, op_2::OneScalar{FT}) where FT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: op_2.tau
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(FT) :: tau12

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)
  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,ibnd]
          op_1.ssa[icol,ilay,igpt] = op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(ε_machine(FT),tau12)
          op_1.tau[icol,ilay,igpt] = tau12
          # p is unchanged
        end
      end
    end
  end
end


# ---------------------------------
# increment nstream by 2stream
function increment_bybnd!(op_1::NStream{FT}, op_2::TwoStream{FT}) where FT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom1, nbnd
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(FT), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: op_1.p
#    real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: op_2.tau, op_2.ssa, op_2.g
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(FT) :: tau12, tauscat12
#    real(FT), dimension(nmom1) :: temp_moms # TK
#    integer  :: imom  #TK

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)
  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,ibnd]
          tauscat12 =
             op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
             op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd]
          #
          # Here assume Henyey-Greenstein
          #
          temp_moms[1] = op_2.g[icol,ilay,ibnd]
          for imom = 2:nmom1
            temp_moms[imom] = temp_moms[imom-1] * op_2.g[icol,ilay,ibnd]
          end
          op_1.p[1:nmom1, icol,ilay,igpt] =
              (op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.p[1:nmom1, icol,ilay,igpt] +
               op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd] * temp_moms[1:nmom1]  ) / max(ε_machine(FT),tauscat12)
          op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(FT),tau12)
          op_1.tau[icol,ilay,igpt] = tau12
        end
      end
    end
  end
end


# ---------------------------------
# increment nstream by nstream
function increment_bybnd!(op_1::NStream{FT}, op_2::NStream{FT}) where FT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom1, nmom2, nbnd
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(FT), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: op_1.p
#    real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: op_2.tau, op_2.ssa
#    real(FT), dimension(nmom2,
#                        ncol,nlay,nbnd), intent(in   ) :: op_2.p
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd, mom_lim
#    real(FT) :: tau12, tauscat12

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)
  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  mom_lim = min(nmom1, nmom2)
  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,ibnd]
          tauscat12 =
             op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
             op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd]
          #
          # If op2 has more moments than op1 these are ignored;
          #   if it has fewer moments the higher orders are assumed to be 0
          #
          op_1.p[1:mom_lim, icol,ilay,igpt] =
              (op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.p[1:mom_lim, icol,ilay,igpt] +
               op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd] * op_2.p[1:mom_lim, icol,ilay,ibnd]) / max( ε_machine(FT),tauscat12)
          op_1.ssa[icol,ilay,igpt] = tauscat12 / max(eps,tau12)
          op_1.tau[icol,ilay,igpt] = tau12
        end
      end
    end
  end
end
