export DensityMatrix, density_matrix, ρ

"""
    DensityMatrix{B, T, MT}

Density Matrix.

- `B`: batch size
- `T`: element type
"""
struct DensityMatrix{B, T, MT<:AbstractArray{T, 3}} <: AbstractRegister{B}
    state::MT
end

"""
    DensityMatrix(state::AbstractArray{T, 3})
    DensityMatrix(state::AbstractMatrix{T})

Create a `DensityMatrix` with a state represented by array.
"""
DensityMatrix(state::MT) where {T, MT<:AbstractArray{T, 3}} = DensityMatrix{size(state, 3), T, MT}(state)
DensityMatrix(state::AbstractMatrix) = DensityMatrix(reshape(state, size(state)..., 1))

"""
    state(ρ::DensityMatrix)

Return the raw state of density matrix `ρ`.
"""
state(ρ::DensityMatrix) = ρ.state

YaoBase.nqubits(ρ::DensityMatrix) = log2dim1(state(ρ))
YaoBase.nactive(ρ::DensityMatrix) = nqubits(ρ)
YaoBase.nbatch(dm::DensityMatrix{B}) where B = B
YaoBase.density_matrix(reg::ArrayReg{1}) = DensityMatrix(reg.state * reg.state')
function YaoBase.density_matrix(reg::ArrayReg{B}) where B
    M = size(reg.state, 1)
    s = reshape(reg |> state, M, :, B)
    out = similar(s, M, M, B)
    for b in 1:B
        @inbounds @views out[:,:,b] = s[:,:,b]*s[:,:,b]'
    end
    return DensityMatrix(out)
end

YaoBase.tracedist(dm1::DensityMatrix{B}, dm2::DensityMatrix{B}) where B =
    map(b->trnorm(dm1.state[:,:,b] - dm2.state[:,:,b]), 1:B)

# TODO: use batch_broadcast in the future
"""
    probs(ρ)

Returns the probability distribution from a density matrix `ρ`.
"""
function YaoBase.probs(m::DensityMatrix{B, T}) where {B, T}
    res = zeros(T, size(m.state, 1), B)
    for i in 1:B
        @inbounds res[:,B] = diag(view(m.state, :,:,i))
    end
    return res
end

YaoBase.probs(m::DensityMatrix{1})= diag(view(m.state, :,:,1))
