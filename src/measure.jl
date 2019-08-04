using StatsBase, StaticArrays, BitBasis, Random
export measure,
    measure!,
    measure_remove!,
    measure_collapseto!,
    select,
    select!

function _measure(rng::AbstractRNG, pl::AbstractVector, nshots::Int)
    N = log2i(length(pl))
    sample(rng, basis(BitStr{N}), Weights(pl), nshots)
end

function _measure(rng::AbstractRNG, pl::AbstractMatrix, nshots::Int)
    B = size(pl, 2)
    res = Matrix{BitStr{log2i(length(pl))}}(undef, nshots, B)
    for ib=1:B
        @inbounds res[:,ib] = _measure(rng, view(pl,:,ib), nshots)
    end
    return res
end

YaoBase.measure(rng::AbstractRNG, ::ComputationalBasis, reg::ArrayReg{1}, ::AllLocs; nshots::Int=1) = _measure(rng, reg |> probs, nshots)

function YaoBase.measure(rng::AbstractRNG, ::ComputationalBasis, reg::ArrayReg{B}, ::AllLocs; nshots::Int=1) where B
    pl = dropdims(sum(reg |> rank3 .|> abs2, dims=2), dims=2)
    return _measure(rng, pl, nshots)
end

function YaoBase.measure_remove!(rng::AbstractRNG, ::ComputationalBasis, reg::ArrayReg{B}, ::AllLocs) where B
    state = reg |> rank3
    nstate = similar(reg.state, 1<<nremain(reg), B)
    pl = dropdims(sum(state .|> abs2, dims=2), dims=2)
    res = Vector{BitStr{nactive(reg)}}(undef, B)
    @inbounds for ib = 1:B
        ires = _measure(rng, view(pl, :, ib), 1)[]
        # notice ires is `BitStr` type, can be use as indices directly.
        nstate[:,ib] = view(state, ires,:,ib)./sqrt(pl[ires, ib])
        res[ib] = ires
    end
    reg.state = reshape(nstate,1,:)
    return res
end

function YaoBase.measure!(rng::AbstractRNG, ::ComputationalBasis, reg::ArrayReg{B}, ::AllLocs) where B
    state = reg |> rank3
    nstate = zero(state)
    res = measure_remove!(rng, reg)
    _nstate = reshape(reg.state, :, B)
    for ib in 1:B
        @inbounds nstate[res[ib], :, ib] .= view(_nstate, :,ib)
    end
    reg.state = reshape(nstate, size(state, 1), :)
    return res
end

function YaoBase.measure_collapseto!(rng::AbstractRNG, ::ComputationalBasis, reg::ArrayReg{B}, ::AllLocs; config::Integer=0) where B
    state = rank3(reg)
    M, N, B1 = size(state)
    nstate = zero(state)
    res = measure_remove!(rng, reg)
    nstate[Int(config)+1, :, :] = reshape(reg.state, :, B)
    reg.state = reshape(nstate, M, N*B)
    return res
end

import YaoBase: select, select!
select(r::ArrayReg{B}, bits::AbstractVector{T}) where {B, T<:BitStr} = ArrayReg{B}(r.state[bits, :])
select(r::ArrayReg{B}, bit::BitStr) where B = select(r, [bit])

function select!(r::ArrayReg, bits::AbstractVector{T}) where T<:BitStr
    r.state = r.state[bits, :]
    return r
end

select!(r::ArrayReg, bit::BitStr) = select!(r, [bit])
