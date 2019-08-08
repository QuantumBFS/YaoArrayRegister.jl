# deprecations
@deprecate select!(r::ArrayReg, bit::Integer) select!(r, BitStr64{nactive(r)}(bit))
@deprecate select(r::ArrayReg, bit::Integer) select(r, BitStr64{nactive(r)}(bit))
@deprecate select!(r::ArrayReg, bits::AbstractVector{<:Integer}) select!(r, reinterpret(BitStr64{nactive(r)}, bits))
@deprecate select(r::ArrayReg, bits::AbstractVector{<:Integer}) select(r, reinterpret(BitStr64{nactive(r)}, bits))
