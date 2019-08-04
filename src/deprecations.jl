# deprecations
@deprecate select!(r::ArrayReg, bit::Integer) select!(r, BitStr{nactive(r)}(bit))
@deprecate select(r::ArrayReg, bit::Integer) select(r, BitStr{nactive(r)}(bit))
@deprecate select!(r::ArrayReg, bits::AbstractVector{<:Integer}) select!(r, reinterpret(BitStr{nactive(r)}, bits))
@deprecate select(r::ArrayReg, bits::AbstractVector{<:Integer}) select(r, reinterpret(BitStr{nactive(r)}, bits))
