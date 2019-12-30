#### Kernels for arrays of optical properties
#  - delta-scaling
#  - adding two sets of properties
#  - extracting subsets
#  - validity checking

ε_machine(FT) = 3*eps(FT)

"""
    delta_scale_kernel!(op::TwoStream{FT}, f) where {FT<:AbstractFloat}

Delta-scale optical properties given

 - `op` optical properties, see [`TwoStream`](@ref)
 - `f` forward scattering [ncol, nlay, ngpt]
"""
function delta_scale_kernel!(op::TwoStream{FT}, f) where {FT<:AbstractFloat}
  @inbounds for igpt = 1:get_ngpt(op)
    @inbounds for ilay = 1:get_nlay(op)
      @inbounds for icol = 1:get_ncol(op)
        wf = op.ssa[icol,ilay,igpt] * f[icol,ilay,igpt]
        op.τ[icol,ilay,igpt] = (FT(1) - wf) * op.τ[icol,ilay,igpt]
        op.ssa[icol,ilay,igpt] = (op.ssa[icol,ilay,igpt] - wf) /  max(ε_machine(FT),(FT(1) - wf))
        op.g[icol,ilay,igpt] = (g[icol,ilay,igpt] - f[icol,ilay,igpt]) /
                                      max(ε_machine(FT),(FT(1) - f[icol,ilay,igpt]))
      end
    end
  end

end

"""
    delta_scale_kernel!(op::TwoStream{FT}) where {FT<:AbstractFloat}

Delta-scale optical properties ``f = g*g`` given

 - `op` optical properties, see [`AbstractOpticalProps`](@ref)
"""
function delta_scale_kernel!(op::TwoStream{FT}) where {FT<:AbstractFloat}
  @inbounds for igpt = 1:get_ngpt(op)
    @inbounds for ilay = 1:get_nlay(op)
      @inbounds for icol = 1:get_ncol(op)
        f  = op.g[icol,ilay,igpt] * op.g[icol,ilay,igpt]
        wf = op.ssa[icol,ilay,igpt] * f
        op.τ[icol,ilay,igpt] = (FT(1) - wf) * op.τ[icol,ilay,igpt]
        op.ssa[icol,ilay,igpt] = (op.ssa[icol,ilay,igpt] - wf) /  max(ε_machine(FT),(FT(1) - wf))
        op.g[icol,ilay,igpt] = (op.g[icol,ilay,igpt] -  f) /  max(ε_machine(FT),(FT(1) -  f))
      end
    end
  end
end

#####
##### Addition of optical properties: the first set are incremented by the second set.
#####

# There are 2 possible representations of optical properties:
#    - `OneScalar` optical depth only
#    - `TwoStream` optical depth, single-scattering albedo, and asymmetry factor g
# Thus we need 4 routines (1 for every combination)

"""
    increment_by_gpoint!(op_1::AbstractOpticalPropsArry{FT},
                         op_2::AbstractOpticalPropsArry{FT}) where {FT<:AbstractFloat}

Adds optical properties by g-point

!!! Note:
    optical properties are defined at the same spectral resolution
"""
increment_by_gpoint!(op_1::OneScalar{FT}, op_2::OneScalar{FT}) where {FT} = (op_1.τ .= op_1.τ .+ op_2.τ)
increment_by_gpoint!(op_1::OneScalar{FT}, op_2::TwoStream{FT}) where {FT} = (op_1.τ .= op_1.τ .+ op_2.τ .* ( FT(1) .- op_2.ssa ))
function increment_by_gpoint!(op_1::TwoStream{FT}, op_2::OneScalar{FT}) where {FT<:AbstractFloat}
  @inbounds for igpt = 1:get_ngpt(op_1)
    @inbounds for ilay = 1:get_nlay(op_1)
      @inbounds for icol = 1:get_ncol(op_1)
        τ12 = op_1.τ[icol,ilay,igpt] + op_2.τ[icol,ilay,igpt]
        op_1.ssa[icol,ilay,igpt] = op_1.τ[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(ε_machine(FT),τ12)
        op_1.τ[icol,ilay,igpt] = τ12
        # g is unchanged
      end
    end
  end
end
function increment_by_gpoint!(op_1::TwoStream{FT}, op_2::TwoStream{FT}) where {FT<:AbstractFloat}
  @inbounds for igpt = 1:get_ngpt(op_1)
    @inbounds for ilay = 1:get_nlay(op_1)
      @inbounds for icol = 1:get_ncol(op_1)
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

"""
    increment_bybnd!(op_1::AbstractOpticalPropsArry{FT},
                     op_2::AbstractOpticalPropsArry{FT}) where {FT<:AbstractFloat}

Adds optical properties by band

!!! Note:
    Assumes second set of optical properties is defined at lower spectral
    resolution (e.g. by band instead of by gpoint)
"""
function increment_bybnd!(op_1::OneScalar{FT}, op_2::OneScalar{FT}) where {FT<:AbstractFloat}
  @inbounds for ibnd = 1:get_nband(op_1)
    @inbounds for igpt = gpt_range(op_1, ibnd)
      op_1.τ[:,:,igpt] .= op_1.τ[:,:,igpt] .+ op_2.τ[:,:,ibnd]
    end
  end
end
function increment_bybnd!(op_1::OneScalar{FT}, op_2::TwoStream{FT}) where {FT<:AbstractFloat}
  @inbounds for ibnd = 1:get_nband(op_1)
    @inbounds for igpt = gpt_range(op_1, ibnd)
      op_1.τ[:,:,igpt] .= op_1.τ[:,:,igpt] .+ op_2.τ[:,:,ibnd] .* ( FT(1) .- op_2.ssa[:,:,ibnd])
    end
  end
end
function increment_bybnd!(op_1::TwoStream{FT}, op_2::OneScalar{FT}) where {FT<:AbstractFloat}
  @inbounds for ibnd = 1:get_nband(op_1)
    @inbounds for igpt = gpt_range(op_1, ibnd)
      @inbounds for ilay = 1:get_nlay(op_1)
        @inbounds for icol = 1:get_ncol(op_1)
          τ12 = op_1.τ[icol,ilay,igpt] + op_2.τ[icol,ilay,ibnd]
          op_1.ssa[icol,ilay,igpt] = op_1.τ[icol,ilay,igpt] * op_1.ssa[icol,ilay,igpt] / max(eps,τ12)
          op_1.τ[icol,ilay,igpt] = τ12
          # g is unchanged
        end
      end
    end
  end
end
function increment_bybnd!(op_1::TwoStream{FT}, op_2::TwoStream{FT}) where {FT<:AbstractFloat}
  @inbounds for ibnd = 1:get_nband(op_1)
    @inbounds for igpt = gpt_range(op_1, ibnd)
      @inbounds for ilay = 1:get_nlay(op_1)
        @inbounds for icol = 1:get_ncol(op_1)
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

