#### Kernels for arrays of optical properties
#  - delta-scaling
#  - adding two sets of properties
#  - extracting subsets
#  - validity checking

ε_machine(FT) = 3*eps(FT)


#
# Delta-scaling, provided only for two-stream properties at present
#

# Delta-scale two-stream optical properties
#   user-provided value of f (forward scattering)
#

"""

#    integer,                               intent(in   ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol, nlay, ngpt), intent(inout) ::  τ, ssa, g
#    real(FT), dimension(ncol, nlay, ngpt), intent(in   ) ::  f

#    real(FT) :: wf
#    integer  :: icol, ilay, igpt
"""
function delta_scale_2str_kernel!(op::AbstractOpticalProps{FT}, f) where FT
  for igpt = 1:get_ngpt(op)
    for ilay = 1:get_nlay(op)
      for icol = 1:get_ncol(op)
        wf = op.ssa[icol,ilay,igpt] * f[icol,ilay,igpt]
        op.τ[icol,ilay,igpt] = (FT(1) - wf) * op.τ[icol,ilay,igpt]
        op.ssa[icol,ilay,igpt] = (op.ssa[icol,ilay,igpt] - wf) /  max(ε_machine(FT),(FT(1) - wf))
        op.g[icol,ilay,igpt] = (g[icol,ilay,igpt] - f[icol,ilay,igpt]) /
                                      max(ε_machine(FT),(FT(1) - f[icol,ilay,igpt]))
      end
    end
  end

end
# ---------------------------------
# Delta-scale
#   f = g*g
#
#    integer,                               intent(in   ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol, nlay, ngpt), intent(inout) ::  τ, ssa, g

#    real(FT) :: f, wf
#    integer  :: icol, ilay, igpt
function delta_scale_2str_kernel!(op::AbstractOpticalProps{FT}) where FT
  for igpt = 1:get_ngpt(op)
    for ilay = 1:get_nlay(op)
      for icol = 1:get_ncol(op)
        f  = op.g[icol,ilay,igpt] * op.g[icol,ilay,igpt]
        wf = op.ssa[icol,ilay,igpt] * f
        op.τ[icol,ilay,igpt] = (FT(1) - wf) * op.τ[icol,ilay,igpt]
        op.ssa[icol,ilay,igpt] = (op.ssa[icol,ilay,igpt] - wf) /  max(ε_machine(FT),(FT(1) - wf))
        op.g[icol,ilay,igpt] = (op.g[icol,ilay,igpt] -  f) /  max(ε_machine(FT),(FT(1) -  f))
      end
    end
  end
end
# -------------------------------------------------------------------------------------------------
#
# Addition of optical properties: the first set are incremented by the second set.
#
#   There are three possible representations of optical properties (scalar = optical depth only;
#   two-stream = τ, single-scattering albedo, and asymmetry factor g, and
#   n-stream = τ, ssa, and phase function moments p.) Thus we need nine routines, three for
#   each choice of representation on the left hand side times three representations of the
#   optical properties to be added.
#
#   There are two sets of these nine routines. In the first the two sets of optical
#   properties are defined at the same spectral resolution. There is also a set of routines
#   to add properties defined at lower spectral resolution to a set defined at higher spectral
#   resolution (adding properties defined by band to those defined by g-point)
#
# -------------------------------------------------------------------------------------------------
#    integer,                              intent(in  ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.τ
#    real(FT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.τ
#    integer  :: icol, ilay, igpt
function increment_by_gpoint!(op_1::OneScalar{FT}, op_2::OneScalar{FT}) where FT
  op_1.τ .= op_1.τ .+ op_2.τ
end
# ---------------------------------
# increment 1scalar by 2stream
function increment_by_gpoint!(op_1::OneScalar{FT}, op_2::TwoStream{FT}) where FT
  op_1.τ .= op_1.τ .+ op_2.τ .* ( FT(1) .- op_2.ssa )
end
# ---------------------------------
# ---------------------------------
# increment 2stream by 1scalar
function increment_by_gpoint!(op_1::TwoStream{FT}, op_2::OneScalar{FT}) where FT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.τ, op_1.ssa
#    real(FT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.τ

#    integer  :: icol, ilay, igpt
#    real(FT) :: τ12

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)

  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = 1:ncol
        τ12 = op_1.τ[icol,ilay,igpt] + op_2.τ[icol,ilay,igpt]
        op_1.ssa[icol,ilay,igpt] = op_1.τ[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(ε_machine(FT),τ12)
        op_1.τ[icol,ilay,igpt] = τ12
        # g is unchanged
      end
    end
  end
end
# ---------------------------------
# increment 2stream by 2stream
function increment_by_gpoint!(op_1::TwoStream{FT}, op_2::TwoStream{FT}) where FT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.τ, op_1.ssa, op_1.g
#    real(FT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.τ, op_2.ssa, op_2.g

#    integer :: icol, ilay, igpt
#    real(FT) :: τ12, τscat12

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)

  for igpt = 1:ngpt
    for ilay = 1:nlay
      for icol = 1:ncol
        # t=op_1.τ + op_2.τ
        τ12 = op_1.τ[icol,ilay,igpt] + op_2.τ[icol,ilay,igpt]
        # w=(op_1.τ*op_1.ssa + op_2.τ*op_2.ssa) / t
        τscat12 = op_1.τ[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
                    op_2.τ[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt]
        op_1.g[icol,ilay,igpt] =
          (op_1.τ[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.g[icol,ilay,igpt] +
           op_2.τ[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt] * op_2.g[icol,ilay,igpt]) / max(ε_machine(FT),τscat12)
        op_1.ssa[icol,ilay,igpt] = τscat12 / max(ε_machine(FT),τ12)
        op_1.τ[icol,ilay,igpt] = τ12
      end
    end
  end
end

#
# Incrementing when the second set of optical properties is defined at lower spectral resolution
#   (e.g. by band instead of by gpoint)
#
function increment_bybnd!(op_1::OneScalar{FT}, op_2::OneScalar{FT}) where FT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: τ1
#    real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: τ2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      op_1.τ[:,:,igpt] .= op_1.τ[:,:,igpt] .+ op_2.τ[:,:,ibnd]
    end
  end
end
# ---------------------------------
# increment 1scalar by 2stream
function increment_bybnd!(op_1::OneScalar{FT}, op_2::TwoStream{FT}) where FT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: τ1
#    real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: τ2, ssa2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)
  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      op_1.τ[:,:,igpt] .= op_1.τ[:,:,igpt] .+ op_2.τ[:,:,ibnd] .* ( FT(1) .- op_2.ssa[:,:,ibnd])
    end
  end
end

"""
    increment_bybnd!(op_1::TwoStream{FT}, op_2::OneScalar{FT})

increment 2stream by 1scalar
integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: τ1, ssa1
real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: τ2
integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

integer  :: icol, ilay, igpt, ibnd
real(FT) :: τ12
"""
function increment_bybnd!(op_1::TwoStream{FT}, op_2::OneScalar{FT}) where FT
  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)
  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      for ilay = 1:nlay
        for icol = 1:ncol
          τ12 = op_1.τ[icol,ilay,igpt] + op_2.τ[icol,ilay,ibnd]
          op_1.ssa[icol,ilay,igpt] = op_1.τ[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(eps,τ12)
          op_1.τ[icol,ilay,igpt] = τ12
          # g is unchanged
        end
      end
    end
  end
end

"""
    increment_bybnd!(op_1::TwoStream{FT}, op_2::TwoStream{FT}) where FT

increment 2stream by 2stream
integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
real(FT), dimension(ncol,nlay,ngpt), intent(inout) :: τ1, ssa1, op_1.g
real(FT), dimension(ncol,nlay,nbnd), intent(in   ) :: τ2, ssa2, op_2.g
integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

integer  :: icol, ilay, igpt, ibnd
real(FT) :: τ12, τscat12
"""
function increment_bybnd!(op_1::TwoStream{FT}, op_2::TwoStream{FT}) where FT

  ngpt = get_ngpt(op_1)
  nlay = get_nlay(op_1)
  ncol = get_ncol(op_1)
  gpt_lims = get_band_lims_gpoint(op_1)
  nbnd = get_nband(op_1)

  for ibnd = 1:nbnd
    for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
      for ilay = 1:nlay
        for icol = 1:ncol
          # t=τ1 + τ2
          τ12 = op_1.τ[icol,ilay,igpt] + op_2.τ[icol,ilay,ibnd]
          # w=(τ1*ssa1 + τ2*ssa2) / t
          τscat12 =
             op_1.τ[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
             op_2.τ[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd]
          op_1.g[icol,ilay,igpt] =
            (op_1.τ[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.g[icol,ilay,igpt] +
             op_2.τ[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd] * op_2.g[icol,ilay,ibnd]) / max(ε_machine(FT),τscat12)
          op_1.ssa[icol,ilay,igpt] = τscat12 / max(ε_machine(FT),τ12)
          op_1.τ[icol,ilay,igpt] = τ12
        end
      end
    end
  end
end

