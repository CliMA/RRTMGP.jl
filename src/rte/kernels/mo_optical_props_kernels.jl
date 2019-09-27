# This code is part of Radiative Transfer for Energetics (RTE)
#
# Contacts: Robert Pincus and Eli Mlawer
# email:  rrtmgp@aer.com
#
# Copyright 2015-2018,  Atmospheric and Environmental Research and
# Regents of the University of Colorado.  All right reserved.
#
# Use and duplication is permitted under the terms of the
#    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
# -------------------------------------------------------------------------------------------------
#
# Kernels for arrays of optical properties:
#   delta-scaling
#   adding two sets of properties
#   extracting subsets
#   validity checking
#
# -------------------------------------------------------------------------------------------------

# module mo_optical_props_kernels

# export delta_scale_2str_f_k,
#        delta_scale_2str_k,
#        increment_1scalar_by_1scalar,
#        increment_1scalar_by_2stream,
#        increment_1scalar_by_nstream,
#        increment_2stream_by_1scalar,
#        increment_2stream_by_2stream,
#        increment_2stream_by_nstream,
#        increment_nstream_by_1scalar,
#        increment_nstream_by_2stream,
#        increment_nstream_by_nstream,
#        inc_1scalar_by_1scalar_bybnd,
#        inc_1scalar_by_2stream_bybnd,
#        inc_1scalar_by_nstream_bybnd,
#        inc_2stream_by_1scalar_bybnd,
#        inc_2stream_by_2stream_bybnd,
#        inc_2stream_by_nstream_bybnd,
#        inc_nstream_by_1scalar_bybnd,
#        inc_nstream_by_2stream_bybnd,
#        extract_subset_dim1_3d,
#        extract_subset_dim2_4d,
#        extract_subset_absorption_tau

ε_machine(DT) = 3*eps(DT)

#  use, intrinsic :: iso_c_binding
#  use mo_rte_kind, only: DT, wl
#  implicit none

#  public
#  interface delta_scale_2str_kernel
#    module procedure delta_scale_2str_f_k, delta_scale_2str_k
#  end interface

#  interface extract_subset
#    module procedure extract_subset_dim1_3d, extract_subset_dim2_4d
#    module procedure extract_subset_absorption_tau
#  end interface extract_subset

#  real(DT), parameter, private :: eps = 3.0_DT*tiny(1.0_DT)
#contains
  # -------------------------------------------------------------------------------------------------
  #
  # Delta-scaling, provided only for two-stream properties at present
  #
  # -------------------------------------------------------------------------------------------------
  # Delta-scale two-stream optical properties
  #   user-provided value of f (forward scattering)
  #
  function delta_scale_2str_kernel!(op::ty_optical_props{DT}, f) where DT
#      bind(C, name="delta_scale_2str_f_k")
#    integer,                               intent(in   ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol, nlay, ngpt), intent(inout) ::  tau, ssa, g
#    real(DT), dimension(ncol, nlay, ngpt), intent(in   ) ::  f

#    real(DT) :: wf
#    integer  :: icol, ilay, igpt

    for igpt = 1:get_ngpt(op)
      for ilay = 1:get_nlay(op)
        for icol = 1:get_ncol(op)
          wf = op.ssa[icol,ilay,igpt] * f[icol,ilay,igpt]
          op.tau[icol,ilay,igpt] = (DT(1) - wf) * op.tau[icol,ilay,igpt]
          op.ssa[icol,ilay,igpt] = (op.ssa[icol,ilay,igpt] - wf) /  max(ε_machine(DT),(DT(1) - wf))
          op.g[icol,ilay,igpt] = (g[icol,ilay,igpt] - f[icol,ilay,igpt]) /
                                        max(ε_machine(DT),(DT(1) - f[icol,ilay,igpt]))
        end
      end
    end

  end
  # ---------------------------------
  # Delta-scale
  #   f = g*g
  #
  function delta_scale_2str_kernel!(op::ty_optical_props{DT}) where DT
#      bind(C, name="delta_scale_2str_k")
#    integer,                               intent(in   ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol, nlay, ngpt), intent(inout) ::  tau, ssa, g

#    real(DT) :: f, wf
#    integer  :: icol, ilay, igpt

    for igpt = 1:get_ngpt(op)
      for ilay = 1:get_nlay(op)
        for icol = 1:get_ncol(op)
          f  = op.g[icol,ilay,igpt] * op.g[icol,ilay,igpt]
          wf = op.ssa[icol,ilay,igpt] * f
          op.tau[icol,ilay,igpt] = (DT(1) - wf) * op.tau[icol,ilay,igpt]
          op.ssa[icol,ilay,igpt] = (op.ssa[icol,ilay,igpt] - wf) /  max(ε_machine(DT),(DT(1) - wf))
          op.g[icol,ilay,igpt] = (op.g[icol,ilay,igpt] -  f) /  max(ε_machine(DT),(DT(1) -  f))
        end
      end
    end

  end
  # -------------------------------------------------------------------------------------------------
  #
  # Addition of optical properties: the first set are incremented by the second set.
  #
  #   There are three possible representations of optical properties (scalar = optical depth only;
  #   two-stream = tau, single-scattering albedo, and asymmetry factor g, and
  #   n-stream = tau, ssa, and phase function moments p.) Thus we need nine routines, three for
  #   each choice of representation on the left hand side times three representations of the
  #   optical properties to be added.
  #
  #   There are two sets of these nine routines. In the first the two sets of optical
  #   properties are defined at the same spectral resolution. There is also a set of routines
  #   to add properties defined at lower spectral resolution to a set defined at higher spectral
  #   resolution (adding properties defined by band to those defined by g-point)
  #
  # -------------------------------------------------------------------------------------------------
  function increment_by_gpoint!(op_1::ty_optical_props_1scl{DT}, op_2::ty_optical_props_1scl{DT}) where DT
  # bind(C, name="increment_1scalar_by_1scalar")
#    integer,                              intent(in  ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau
#    real(DT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau

#    integer  :: icol, ilay, igpt

    op_1.tau .= op_1.tau .+ op_2.tau
  end
  # ---------------------------------
  # increment 1scalar by 2stream
  function increment_by_gpoint!(op_1::ty_optical_props_1scl{DT}, op_2::ty_optical_props_2str{DT}) where DT
    op_1.tau .= op_1.tau .+ op_2.tau .* ( DT(1) .- op_2.ssa[icol,ilay,igpt])
  end
  # ---------------------------------
  # increment 1scalar by nstream
  function increment_by_gpoint!(op_1::ty_optical_props_1scl{DT}, op_2::ty_optical_props_nstr{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau
#    real(DT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa

#    integer  :: icol, ilay, igpt

    op_1.tau = op_1.tau .+ op_2.tau * ( DT(1) - op_2.ssa)
  end
  # ---------------------------------
  # ---------------------------------
  # increment 2stream by 1scalar
  function increment_by_gpoint!(op_1::ty_optical_props_2str{DT}, op_2::ty_optical_props_1scl{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(DT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau

#    integer  :: icol, ilay, igpt
#    real(DT) :: tau12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,igpt]
          op_1.ssa[icol,ilay,igpt] = op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(ε_machine(DT),tau12)
          op_1.tau[icol,ilay,igpt] = tau12
          # g is unchanged
        end
      end
    end
  end
  # ---------------------------------
  # increment 2stream by 2stream
  function increment_by_gpoint!(op_1::ty_optical_props_2str{DT}, op_2::ty_optical_props_2str{DT}) where DT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa, op_1.g
#    real(DT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa, op_2.g

#    integer :: icol, ilay, igpt
#    real(DT) :: tau12, tauscat12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          # t=op_1.tau + op_2.tau
          tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,igpt]
          # w=(op_1.tau*op_1.ssa + op_2.tau*op_2.ssa) / t
          tauscat12 = op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] +
                      op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt]
          op_1.g[icol,ilay,igpt] =
            (op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] * op_1.g[icol,ilay,igpt] +
             op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt] * op_2.g[icol,ilay,igpt]) / max(ε_machine(DT),tauscat12)
          op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(DT),tau12)
          op_1.tau[icol,ilay,igpt] = tau12
        end
      end
    end
  end
  # ---------------------------------
  # increment 2stream by nstream
  function increment_by_gpoint!(op_1::ty_optical_props_2str{DT}, op_2::ty_optical_props_nstr{DT}) where DT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom2
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa, op_1.g
#    real(DT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa
#    real(DT), dimension(nmom2,
#                        ncol,nlay,ngpt), intent(in   ) :: op_2.p

#    integer  :: icol, ilay, igpt
#    real(DT) :: tau12, tauscat12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)

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
             op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt] * op_2.p[1, icol,ilay,igpt]) / max(ε_machine(DT),tauscat12)
          op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(DT),tau12)
          op_1.tau[icol,ilay,igpt] = tau12
        end
      end
    end
  end
  # ---------------------------------
  # ---------------------------------
  # increment nstream by 1scalar
  function increment_by_gpoint!(op_1::ty_optical_props_nstr{DT}, op_2::ty_optical_props_1scl{DT}) where DT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(DT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau

#    integer  :: icol, ilay, igpt
#    real(DT) :: tau12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,igpt]
          op_1.ssa[icol,ilay,igpt] = op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(ε_machine(DT),tau12)
          op_1.tau[icol,ilay,igpt] = tau12
          # p is unchanged
        end
      end
    end
  end
  # ---------------------------------
  # increment nstream by 2stream
  function increment_by_gpoint!(op_1::ty_optical_props_nstr{DT}, op_2::ty_optical_props_2str{DT}) where DT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom1
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(DT), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: op_1.p
#    real(DT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa, op_2.g

#    integer  :: icol, ilay, igpt
#    real(DT) :: tau12, tauscat12
#    real(DT), dimension(nmom1) :: temp_moms # TK
#    integer  :: imom  #TK

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)

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
               op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt] * temp_moms[1:nmom1]  ) / max(ε_machine(DT),tauscat12)
          op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(DT),tau12)
          op_1.tau[icol,ilay,igpt] = tau12
        end
      end
    end
  end
  # ---------------------------------
  # increment nstream by nstream
  function increment_by_gpoint!(op_1::ty_optical_props_nstr{DT}, op_2::ty_optical_props_nstr{DT}) where DT
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom1, nmom2
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(DT), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: op_1.p
#    real(DT), dimension(ncol,nlay,ngpt), intent(in   ) :: op_2.tau, op_2.ssa
#    real(DT), dimension(nmom2,
#                        ncol,nlay,ngpt), intent(in   ) :: op_2.p

#    integer  :: icol, ilay, igpt, mom_lim
#    real(DT) :: tau12, tauscat12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)

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
               op_2.tau[icol,ilay,igpt] * op_2.ssa[icol,ilay,igpt] * op_2.p[1:mom_lim, icol,ilay,igpt]) / max( ε_machine(DT),tauscat12)
          op_1.ssa[icol,ilay,igpt] = tauscat12 / max(eps,tau12)
          op_1.tau[icol,ilay,igpt] = tau12
        end
      end
    end
  end
  # -------------------------------------------------------------------------------------------------
  #
  # Incrementing when the second set of optical properties is defined at lower spectral resolution
  #   (e.g. by band instead of by gpoint)
  #
  # -------------------------------------------------------------------------------------------------
  function increment_bybnd!(op_1::ty_optical_props_1scl{DT}, op_2::ty_optical_props_1scl{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        op_1.tau[:,:,igpt] .= op_1.tau[:,:,igpt] .+ op_2.tau[:,:,ibnd]
      end
    end
  end
  # ---------------------------------
  # increment 1scalar by 2stream
  function increment_bybnd!(op_1::ty_optical_props_1scl{DT}, op_2::ty_optical_props_2str{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)
    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        op_1.tau[:,:,igpt] .= op_1.tau[:,:,igpt] .+ op_2.tau[:,:,ibnd] .* ( DT(1) .- op_2.ssa[:,:,ibnd])
      end
    end
  end
  # ---------------------------------
  # increment 1scalar by nstream
  function increment_bybnd!(op_1::ty_optical_props_1scl{DT}, op_2::ty_optical_props_nstr{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)
    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        op_1.tau[:,:,igpt] .= op_1.tau[:,:,igpt] .+ op_2.tau[:,:,ibnd] .* ( DT(1) .- op_2.ssa[:,:,ibnd])
      end
    end
  end

    # ---------------------------------
  # increment 2stream by 1scalar
  function increment_bybnd!(op_1::ty_optical_props_2str{DT}, op_2::ty_optical_props_1scl{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(DT) :: tau12
    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)
    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        for ilay = 1:nlay
          for icol = 1:ncol
            tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,ibnd]
            op_1.ssa[icol,ilay,igpt] = op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(eps,tau12)
            op_1.tau[icol,ilay,igpt] = tau12
            # g is unchanged
          end
        end
      end
    end
  end
  # ---------------------------------
  # increment 2stream by 2stream
  function increment_bybnd!(op_1::ty_optical_props_2str{DT}, op_2::ty_optical_props_2str{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, op_1.g
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2, op_2.g
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(DT) :: tau12, tauscat12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)
    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        for ilay = 1:nlay
          for icol = 1:ncol
            # t=tau1 + tau2
            tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,ibnd]
            # w=(tau1*ssa1 + tau2*ssa2) / t
            tauscat12 =
               tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] +
               tau2[icol,ilay,ibnd] * ssa2[icol,ilay,ibnd]
            op_1.g[icol,ilay,igpt] =
              (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * op_1.g[icol,ilay,igpt] +
               tau2[icol,ilay,ibnd] * ssa2[icol,ilay,ibnd] * op_2.g[icol,ilay,ibnd]) / max(ε_machine(DT),tauscat12)
            ssa1[icol,ilay,igpt] = tauscat12 / max(ε_machine(DT),tau12)
            tau1[icol,ilay,igpt] = tau12
          end
        end
      end
    end
  end
  # ---------------------------------
  # increment 2stream by nstream
  function increment_bybnd!(op_1::ty_optical_props_2str{DT}, op_2::ty_optical_props_nstr{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom2, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa, op_1.g
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: op_2.tau, op_2.ssa
#    real(DT), dimension(nmom2,
#                        ncol,nlay,nbnd), intent(in   ) :: op_2.p
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(DT) :: tau12, tauscat12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)
    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

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
               op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd] * op_2.p[1, icol,ilay,ibnd]) / max(ε_machine(DT),tauscat12)
            op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(DT),tau12)
            op_1.tau[icol,ilay,igpt] = tau12
          end
        end
      end
    end
  end
  # ---------------------------------
  # ---------------------------------
  # increment nstream by 1scalar
  function increment_bybnd!(op_1::ty_optical_props_nstr{DT}, op_2::ty_optical_props_1scl{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: op_2.tau
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(DT) :: tau12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)
    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        for ilay = 1:nlay
          for icol = 1:ncol
            tau12 = op_1.tau[icol,ilay,igpt] + op_2.tau[icol,ilay,ibnd]
            op_1.ssa[icol,ilay,igpt] = op_1.tau[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(ε_machine(DT),tau12)
            op_1.tau[icol,ilay,igpt] = tau12
            # p is unchanged
          end
        end
      end
    end
  end
  # ---------------------------------
  # increment nstream by 2stream
  function increment_bybnd!(op_1::ty_optical_props_nstr{DT}, op_2::ty_optical_props_2str{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom1, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(DT), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: op_1.p
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: op_2.tau, op_2.ssa, op_2.g
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(DT) :: tau12, tauscat12
#    real(DT), dimension(nmom1) :: temp_moms # TK
#    integer  :: imom  #TK

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)
    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

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
                 op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd] * temp_moms[1:nmom1]  ) / max(ε_machine(DT),tauscat12)
            op_1.ssa[icol,ilay,igpt] = tauscat12 / max(ε_machine(DT),tau12)
            op_1.tau[icol,ilay,igpt] = tau12
          end
        end
      end
    end
  end
  # ---------------------------------
  # increment nstream by nstream
  function increment_bybnd!(op_1::ty_optical_props_nstr{DT}, op_2::ty_optical_props_nstr{DT}) where DT
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom1, nmom2, nbnd
#    real(DT), dimension(ncol,nlay,ngpt), intent(inout) :: op_1.tau, op_1.ssa
#    real(DT), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: op_1.p
#    real(DT), dimension(ncol,nlay,nbnd), intent(in   ) :: op_2.tau, op_2.ssa
#    real(DT), dimension(nmom2,
#                        ncol,nlay,nbnd), intent(in   ) :: op_2.p
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd, mom_lim
#    real(DT) :: tau12, tauscat12

    ngpt = get_ngpt(op_io)
    nlay = get_nlay(op_io)
    ncol = get_ncol(op_io)
    gpt_lims = get_band_lims_gpoint(op_1)
    nbnd = get_nband(op_io)

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
                 op_2.tau[icol,ilay,ibnd] * op_2.ssa[icol,ilay,ibnd] * op_2.p[1:mom_lim, icol,ilay,ibnd]) / max( ε_machine(DT),tauscat12)
            op_1.ssa[icol,ilay,igpt] = tauscat12 / max(eps,tau12)
            op_1.tau[icol,ilay,igpt] = tau12
          end
        end
      end
    end
  end
  # -------------------------------------------------------------------------------------------------
  #
  # Subsetting, meaning extracting some portion of the 3D domain
  #
  # -------------------------------------------------------------------------------------------------
  function extract_subset_dim1_3d(ncol, nlay, ngpt, array_in, colS, colE, array_out)#
#    bind (C, name="extract_subset_dim1_3d")
#    integer,                             intent(in ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol,nlay,ngpt), intent(in ) :: array_in
#    integer,                             intent(in ) :: colS, colE
#    real(DT), dimension(colE-colS+1,
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
#    real(DT), dimension(nmom,ncol,nlay,ngpt), intent(in ) :: array_in
#    integer,                                  intent(in ) :: colS, colE
#    real(DT), dimension(nmom,colE-colS+1,
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
  # ---------------------------------
  #
  # Extract the absorption optical thickness which requires mulitplying by 1 - ssa
  #
  function extract_subset_absorption_tau(ncol, nlay, ngpt, tau_in, ssa_in, colS, colE, tau_out) # bind (C, name="extract_subset_absorption_tau")
#    integer,                             intent(in ) :: ncol, nlay, ngpt
#    real(DT), dimension(ncol,nlay,ngpt), intent(in ) :: tau_in, ssa_in
#    integer,                             intent(in ) :: colS, colE
#    real(DT), dimension(colE-colS+1,
#                             nlay,ngpt), intent(out) :: tau_out

#    integer :: icol, ilay, igpt

    DT = eltype(tau_in)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = colS:colE
          tau_out[icol-colS+1, ilay, igpt] = tau_in[icol, ilay, igpt] * ( DT(1) - ssa_in[icol, ilay, igpt] )
        end
      end
    end

  end
  #-------------------------------------
# end
