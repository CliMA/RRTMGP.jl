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

module mo_optical_props_kernels

export delta_scale_2str_f_k,
       delta_scale_2str_k,
       increment_1scalar_by_1scalar,
       increment_1scalar_by_2stream,
       increment_1scalar_by_nstream,
       increment_2stream_by_1scalar,
       increment_2stream_by_2stream,
       increment_2stream_by_nstream,
       increment_nstream_by_1scalar,
       increment_nstream_by_2stream,
       increment_nstream_by_nstream,
       inc_1scalar_by_1scalar_bybnd,
       inc_1scalar_by_2stream_bybnd,
       inc_1scalar_by_nstream_bybnd,
       inc_2stream_by_1scalar_bybnd,
       inc_2stream_by_2stream_bybnd,
       inc_2stream_by_nstream_bybnd,
       inc_nstream_by_1scalar_bybnd,
       inc_nstream_by_2stream_bybnd,
       extract_subset_dim1_3d,
       extract_subset_dim2_4d,
       extract_subset_absorption_tau

ε_machine(DT) = 3*eps(DT)

#  use, intrinsic :: iso_c_binding
#  use mo_rte_kind, only: wp, wl
#  implicit none

#  public
#  interface delta_scale_2str_kernel
#    module procedure delta_scale_2str_f_k, delta_scale_2str_k
#  end interface

#  interface extract_subset
#    module procedure extract_subset_dim1_3d, extract_subset_dim2_4d
#    module procedure extract_subset_absorption_tau
#  end interface extract_subset

#  real(wp), parameter, private :: eps = 3.0_wp*tiny(1.0_wp)
#contains
  # -------------------------------------------------------------------------------------------------
  #
  # Delta-scaling, provided only for two-stream properties at present
  #
  # -------------------------------------------------------------------------------------------------
  # Delta-scale two-stream optical properties
  #   user-provided value of f (forward scattering)
  #
  function delta_scale_2str_f_k(ncol, nlay, ngpt, tau, ssa, g, f)
#      bind(C, name="delta_scale_2str_f_k")
#    integer,                               intent(in   ) :: ncol, nlay, ngpt
#    real(wp), dimension(ncol, nlay, ngpt), intent(inout) ::  tau, ssa, g
#    real(wp), dimension(ncol, nlay, ngpt), intent(in   ) ::  f

#    real(wp) :: wf
#    integer  :: icol, ilay, igpt

    wp = eltype(tau)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          wf = ssa[icol,ilay,igpt] * f[icol,ilay,igpt]
          tau[icol,ilay,igpt] = (wp(1) - wf) * tau[icol,ilay,igpt]
          ssa[icol,ilay,igpt] = (ssa[icol,ilay,igpt] - wf) /  max(ε_machine(wp),(wp(1) - wf))
          g[icol,ilay,igpt] = (g[icol,ilay,igpt] - f[icol,ilay,igpt]) /
                                        max(ε_machine(wp),(wp(1) - f[icol,ilay,igpt]))
        end
      end
    end

  end
  # ---------------------------------
  # Delta-scale
  #   f = g*g
  #
  function delta_scale_2str_k(ncol, nlay, ngpt, tau, ssa, g)
#      bind(C, name="delta_scale_2str_k")
#    integer,                               intent(in   ) :: ncol, nlay, ngpt
#    real(wp), dimension(ncol, nlay, ngpt), intent(inout) ::  tau, ssa, g

#    real(wp) :: f, wf
#    integer  :: icol, ilay, igpt

    wp = eltype(tau)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          f  = g[icol,ilay,igpt] * g[icol,ilay,igpt]
          wf = ssa[icol,ilay,igpt] * f
          tau[icol,ilay,igpt] = (wp(1) - wf) * tau[icol,ilay,igpt]
          ssa[icol,ilay,igpt] = (ssa[icol,ilay,igpt] - wf) /  max(ε_machine(wp),(wp(1) - wf))
          g[icol,ilay,igpt] = (g[icol,ilay,igpt] -  f) /  max(ε_machine(wp),(wp(1) -  f))
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
  function increment_1scalar_by_1scalar(ncol, nlay, ngpt,
                                               tau1,
                                               tau2) # bind(C, name="increment_1scalar_by_1scalar")
#    integer,                              intent(in  ) :: ncol, nlay, ngpt
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2

#    integer  :: icol, ilay, igpt

    wp = eltype(tau1)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau1[icol,ilay,igpt] = tau1[icol,ilay,igpt] + tau2[icol,ilay,igpt]
        end
      end
    end
  end
  # ---------------------------------
  # increment 1scalar by 2stream
  function increment_1scalar_by_2stream(ncol, nlay, ngpt,
                                               tau1,
                                               tau2, ssa2) # bind(C, name="increment_1scalar_by_2stream")
#    integer,                              intent(in   ) :: ncol, nlay, ngpt
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2

#    integer  :: icol, ilay, igpt

    wp = eltype(tau1)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau1[icol,ilay,igpt] = tau1[icol,ilay,igpt] +
                                 tau2[icol,ilay,igpt] * ( wp(1) - ssa2[icol,ilay,igpt])
        end
      end
    end
  end
  # ---------------------------------
  # increment 1scalar by nstream
  function increment_1scalar_by_nstream(ncol, nlay, ngpt,
                                               tau1,
                                               tau2, ssa2) # bind(C, name="increment_1scalar_by_nstream")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2

#    integer  :: icol, ilay, igpt

    wp = eltype(tau1)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau1[icol,ilay,igpt] = tau1[icol,ilay,igpt] +
                                 tau2[icol,ilay,igpt] * ( wp(1) - ssa2[icol,ilay,igpt])
        end
      end
    end
  end
  # ---------------------------------
  # ---------------------------------
  # increment 2stream by 1scalar
  function increment_2stream_by_1scalar(ncol, nlay, ngpt,
                                               tau1, ssa1,
                                               tau2) # bind(C, name="increment_2stream_by_1scalar")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2

#    integer  :: icol, ilay, igpt
#    real(wp) :: tau12

    wp = eltype(tau1)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,igpt]
          ssa1[icol,ilay,igpt] = tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] / max(ε_machine(wp),tau12)
          tau1[icol,ilay,igpt] = tau12
          # g is unchanged
        end
      end
    end
  end
  # ---------------------------------
  # increment 2stream by 2stream
  function increment_2stream_by_2stream(ncol, nlay, ngpt,
                                               tau1, ssa1, g1,
                                               tau2, ssa2, g2) # bind(C, name="increment_2stream_by_2stream")
#    integer,                              intent(in   ) :: ncol, nlay, ngpt
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, g1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2, g2

#    integer :: icol, ilay, igpt
#    real(wp) :: tau12, tauscat12
    wp = eltype(tau12)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          # t=tau1 + tau2
          tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,igpt]
          # w=(tau1*ssa1 + tau2*ssa2) / t
          tauscat12 = tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] +
                      tau2[icol,ilay,igpt] * ssa2[icol,ilay,igpt]
          g1[icol,ilay,igpt] =
            (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * g1[icol,ilay,igpt] +
             tau2[icol,ilay,igpt] * ssa2[icol,ilay,igpt] * g2[icol,ilay,igpt]) / max(ε_machine(wp),tauscat12)
          ssa1[icol,ilay,igpt] = tauscat12 / max(ε_machine(wp),tau12)
          tau1[icol,ilay,igpt] = tau12
        end
      end
    end
  end
  # ---------------------------------
  # increment 2stream by nstream
  function increment_2stream_by_nstream(ncol, nlay, ngpt, nmom2,
                                               tau1, ssa1, g1,
                                               tau2, ssa2, p2) #bind(C, name="increment_2stream_by_nstream")
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom2
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, g1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2
#    real(wp), dimension(nmom2,
#                        ncol,nlay,ngpt), intent(in   ) :: p2

#    integer  :: icol, ilay, igpt
#    real(wp) :: tau12, tauscat12

    wp = eltype(tau12)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          # t=tau1 + tau2
          tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,igpt]
          # w=(tau1*ssa1 + tau2*ssa2) / t
          tauscat12 =
             tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] +
             tau2[icol,ilay,igpt] * ssa2[icol,ilay,igpt]
          g1[icol,ilay,igpt] =
            (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * g1[   icol,ilay,igpt]+
             tau2[icol,ilay,igpt] * ssa2[icol,ilay,igpt] * p2[1, icol,ilay,igpt]) / max(ε_machine(wp),tauscat12)
          ssa1[icol,ilay,igpt] = tauscat12 / max(ε_machine(wp),tau12)
          tau1[icol,ilay,igpt] = tau12
        end
      end
    end
  end
  # ---------------------------------
  # ---------------------------------
  # increment nstream by 1scalar
  function increment_nstream_by_1scalar(ncol, nlay, ngpt,
                                               tau1, ssa1,
                                               tau2) # bind(C, name="increment_nstream_by_1scalar")
#    integer,                              intent(in   ) :: ncol, nlay, ngpt
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2

#    integer  :: icol, ilay, igpt
#    real(wp) :: tau12

    wp = eltype(tau12)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,igpt]
          ssa1[icol,ilay,igpt] = tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] / max(ε_machine(wp),tau12)
          tau1[icol,ilay,igpt] = tau12
          # p is unchanged
        end
      end
    end
  end
  # ---------------------------------
  # increment nstream by 2stream
  function increment_nstream_by_2stream(ncol, nlay, ngpt, nmom1,
                                               tau1, ssa1, p1,
                                               tau2, ssa2, g2) #bind(C, name="increment_nstream_by_2stream")
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom1
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(wp), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: p1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2, g2

#    integer  :: icol, ilay, igpt
#    real(wp) :: tau12, tauscat12
#    real(wp), dimension(nmom1) :: temp_moms # TK
#    integer  :: imom  #TK

    wp = eltype(tau12)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,igpt]
          tauscat12 =
             tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] +
             tau2[icol,ilay,igpt] * ssa2[icol,ilay,igpt]
          #
          # Here assume Henyey-Greenstein
          #
          temp_moms[1] = g2[icol,ilay,igpt]
          for imom = 2:nmom1
            temp_moms[imom] = temp_moms[imom-1] * g2[icol,ilay,igpt]
          end
          p1[1:nmom1, icol,ilay,igpt] =
              (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * p1[1:nmom1, icol,ilay,igpt] +
               tau2[icol,ilay,igpt] * ssa2[icol,ilay,igpt] * temp_moms[1:nmom1]  ) / max(ε_machine(wp),tauscat12)
          ssa1[icol,ilay,igpt] = tauscat12 / max(ε_machine(wp),tau12)
          tau1[icol,ilay,igpt] = tau12
        end
      end
    end
  end
  # ---------------------------------
  # increment nstream by nstream
  function increment_nstream_by_nstream(ncol, nlay, ngpt, nmom1, nmom2,
                                               tau1, ssa1, p1,
                                               tau2, ssa2, p2) # bind(C, name="increment_nstream_by_nstream")
#    integer,                              intent(in   ) :: ncol, nlay, ngpt, nmom1, nmom2
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(wp), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: p1
#    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: tau2, ssa2
#    real(wp), dimension(nmom2,
#                        ncol,nlay,ngpt), intent(in   ) :: p2

#    integer  :: icol, ilay, igpt, mom_lim
#    real(wp) :: tau12, tauscat12


    wp = eltype(tau12)

    mom_lim = min(nmom1, nmom2)
    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = 1:ncol
          tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,igpt]
          tauscat12 =
             tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] +
             tau2[icol,ilay,igpt] * ssa2[icol,ilay,igpt]
          #
          # If op2 has more moments than op1 these are ignored;
          #   if it has fewer moments the higher orders are assumed to be 0
          #
          p1[1:mom_lim, icol,ilay,igpt] =
              (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * p1[1:mom_lim, icol,ilay,igpt] +
               tau2[icol,ilay,igpt] * ssa2[icol,ilay,igpt] * p2[1:mom_lim, icol,ilay,igpt]) / max( ε_machine(wp),tauscat12)
          ssa1[icol,ilay,igpt] = tauscat12 / max(eps,tau12)
          tau1[icol,ilay,igpt] = tau12
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
  function inc_1scalar_by_1scalar_bybnd(ncol, nlay, ngpt,
                                               tau1,
                                               tau2,
                                               nbnd, gpt_lims) # bind(C, name="inc_1scalar_by_1scalar_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        tau1[:,:,igpt] = tau1[:,:,igpt] + tau2[:,:,ibnd]
      end
    end
  end
  # ---------------------------------
  # increment 1scalar by 2stream
  function inc_1scalar_by_2stream_bybnd(ncol, nlay, ngpt,
                                               tau1,
                                               tau2, ssa2,
                                               nbnd, gpt_lims) # bind(C, name="inc_1scalar_by_2stream_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

    wp = eltype(tau12)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        tau1[:,:,igpt] = tau1[:,:,igpt] + tau2[:,:,ibnd] * ( wp(1) - ssa2[:,:,ibnd])
      end
    end
  end
  # ---------------------------------
  # increment 1scalar by nstream
  function inc_1scalar_by_nstream_bybnd(ncol, nlay, ngpt,
                                               tau1,
                                               tau2, ssa2,
                                               nbnd, gpt_lims) # bind(C, name="inc_1scalar_by_nstream_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer :: ibnd, igpt

    wp = eltype(tau12)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        tau1[:,:,igpt] = tau1[:,:,igpt] + tau2[:,:,ibnd] * ( wp(1) - ssa2[:,:,ibnd])
      end
    end
  end

    # ---------------------------------
  # increment 2stream by 1scalar
  function inc_2stream_by_1scalar_bybnd(ncol, nlay, ngpt,
                                               tau1, ssa1,
                                               tau2,
                                               nbnd, gpt_lims) # bind(C, name="inc_2stream_by_1scalar_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(wp) :: tau12

    wp = eltype(tau12)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        for ilay = 1:nlay
          for icol = 1:ncol
            tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,ibnd]
            ssa1[icol,ilay,igpt] = tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] / max(eps,tau12)
            tau1[icol,ilay,igpt] = tau12
            # g is unchanged
          end
        end
      end
    end
  end
  # ---------------------------------
  # increment 2stream by 2stream
  function inc_2stream_by_2stream_bybnd(ncol, nlay, ngpt,
                                               tau1, ssa1, g1,
                                               tau2, ssa2, g2,
                                               nbnd, gpt_lims) # bind(C, name="inc_2stream_by_2stream_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, g1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2, g2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(wp) :: tau12, tauscat12

    wp = eltype(tau12)

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
            g1[icol,ilay,igpt] =
              (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * g1[icol,ilay,igpt] +
               tau2[icol,ilay,ibnd] * ssa2[icol,ilay,ibnd] * g2[icol,ilay,ibnd]) / max(ε_machine(wp),tauscat12)
            ssa1[icol,ilay,igpt] = tauscat12 / max(ε_machine(wp),tau12)
            tau1[icol,ilay,igpt] = tau12
          end
        end
      end
    end
  end
  # ---------------------------------
  # increment 2stream by nstream
  function inc_2stream_by_nstream_bybnd(ncol, nlay, ngpt, nmom2,
                                               tau1, ssa1, g1,
                                               tau2, ssa2, p2,
                                               nbnd, gpt_lims) # bind(C, name="inc_2stream_by_nstream_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom2, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1, g1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
#    real(wp), dimension(nmom2,
#                        ncol,nlay,nbnd), intent(in   ) :: p2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(wp) :: tau12, tauscat12

    wp = eltype(tau12)

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
            g1[icol,ilay,igpt] =
              (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * g1[   icol,ilay,igpt]+
               tau2[icol,ilay,ibnd] * ssa2[icol,ilay,ibnd] * p2[1, icol,ilay,ibnd]) / max(ε_machine(wp),tauscat12)
            ssa1[icol,ilay,igpt] = tauscat12 / max(ε_machine(wp),tau12)
            tau1[icol,ilay,igpt] = tau12
          end
        end
      end
    end
  end
  # ---------------------------------
  # ---------------------------------
  # increment nstream by 1scalar
  function inc_nstream_by_1scalar_bybnd(ncol, nlay, ngpt,
                                               tau1, ssa1,
                                               tau2,
                                               nbnd, gpt_lims) #bind(C, name="inc_nstream_by_1scalar_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(wp) :: tau12

    wp = eltype(tau12)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        for ilay = 1:nlay
          for icol = 1:ncol
            tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,ibnd]
            ssa1[icol,ilay,igpt] = tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] / max(ε_machine(wp),tau12)
            tau1[icol,ilay,igpt] = tau12
            # p is unchanged
          end
        end
      end
    end
  end
  # ---------------------------------
  # increment nstream by 2stream
  function inc_nstream_by_2stream_bybnd(ncol, nlay, ngpt, nmom1,
                                               tau1, ssa1, p1,
                                               tau2, ssa2, g2,
                                               nbnd, gpt_lims)# bind(C, name="inc_nstream_by_2stream_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom1, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(wp), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: p1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2, g2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd
#    real(wp) :: tau12, tauscat12
#    real(wp), dimension(nmom1) :: temp_moms # TK
#    integer  :: imom  #TK

    wp = eltype(tau12)

    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        for ilay = 1:nlay
          for icol = 1:ncol
            tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,ibnd]
            tauscat12 =
               tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] +
               tau2[icol,ilay,ibnd] * ssa2[icol,ilay,ibnd]
            #
            # Here assume Henyey-Greenstein
            #
            temp_moms[1] = g2[icol,ilay,ibnd]
            for imom = 2:nmom1
              temp_moms[imom] = temp_moms[imom-1] * g2[icol,ilay,ibnd]
            end
            p1[1:nmom1, icol,ilay,igpt] =
                (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * p1[1:nmom1, icol,ilay,igpt] +
                 tau2[icol,ilay,ibnd] * ssa2[icol,ilay,ibnd] * temp_moms[1:nmom1]  ) / max(ε_machine(wp),tauscat12)
            ssa1[icol,ilay,igpt] = tauscat12 / max(ε_machine(wp),tau12)
            tau1[icol,ilay,igpt] = tau12
          end
        end
      end
    end
  end
  # ---------------------------------
  # increment nstream by nstream
  function inc_nstream_by_nstream_bybnd(ncol, nlay, ngpt, nmom1, nmom2,
                                               tau1, ssa1, p1,
                                               tau2, ssa2, p2,
                                               nbnd, gpt_lims) #bind(C, name="inc_nstream_by_nstream_bybnd")
#    integer,                             intent(in   ) :: ncol, nlay, ngpt, nmom1, nmom2, nbnd
#    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau1, ssa1
#    real(wp), dimension(nmom1,
#                        ncol,nlay,ngpt), intent(inout) :: p1
#    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: tau2, ssa2
#    real(wp), dimension(nmom2,
#                        ncol,nlay,nbnd), intent(in   ) :: p2
#    integer,  dimension(2,nbnd),         intent(in   ) :: gpt_lims # Starting and ending gpoint for each band

#    integer  :: icol, ilay, igpt, ibnd, mom_lim
#    real(wp) :: tau12, tauscat12

    wp = eltype(tau12)

    mom_lim = min(nmom1, nmom2)
    for ibnd = 1:nbnd
      for igpt = gpt_lims[1, ibnd]:gpt_lims[2, ibnd]
        for ilay = 1:nlay
          for icol = 1:ncol
            tau12 = tau1[icol,ilay,igpt] + tau2[icol,ilay,ibnd]
            tauscat12 =
               tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] +
               tau2[icol,ilay,ibnd] * ssa2[icol,ilay,ibnd]
            #
            # If op2 has more moments than op1 these are ignored;
            #   if it has fewer moments the higher orders are assumed to be 0
            #
            p1[1:mom_lim, icol,ilay,igpt] =
                (tau1[icol,ilay,igpt] * ssa1[icol,ilay,igpt] * p1[1:mom_lim, icol,ilay,igpt] +
                 tau2[icol,ilay,ibnd] * ssa2[icol,ilay,ibnd] * p2[1:mom_lim, icol,ilay,ibnd]) / max( ε_machine(wp),tauscat12)
            ssa1[icol,ilay,igpt] = tauscat12 / max(eps,tau12)
            tau1[icol,ilay,igpt] = tau12
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
#    real(wp), dimension(ncol,nlay,ngpt), intent(in ) :: array_in
#    integer,                             intent(in ) :: colS, colE
#    real(wp), dimension(colE-colS+1,
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
#    real(wp), dimension(nmom,ncol,nlay,ngpt), intent(in ) :: array_in
#    integer,                                  intent(in ) :: colS, colE
#    real(wp), dimension(nmom,colE-colS+1,
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
#    real(wp), dimension(ncol,nlay,ngpt), intent(in ) :: tau_in, ssa_in
#    integer,                             intent(in ) :: colS, colE
#    real(wp), dimension(colE-colS+1,
#                             nlay,ngpt), intent(out) :: tau_out

#    integer :: icol, ilay, igpt

    wp = eltype(tau_in)

    for igpt = 1:ngpt
      for ilay = 1:nlay
        for icol = colS:colE
          tau_out[icol-colS+1, ilay, igpt] = tau_in[icol, ilay, igpt] * ( wp(1) - ssa_in[icol, ilay, igpt] )
        end
      end
    end

  end
  #-------------------------------------
end
