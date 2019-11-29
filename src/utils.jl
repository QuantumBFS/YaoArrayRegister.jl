using LuxurySparse
using TupleTools
@static if hasmethod(TupleTools.diff, Tuple{Tuple{}})
    tuple_diff(args...) = TupleTools.diff(args...)
else
    tuple_diff(v::Tuple{}) = () # similar to diff([])
    tuple_diff(v::Tuple{Any}) = ()
    tuple_diff(v::Tuple) = (v[2] - v[1], tuple_diff(Base.tail(v))...)
end

"""
    sort_unitary(U, locations::NTuple{N, Int}) -> U

Return an sorted unitary operator according to the locations.
"""
function sort_unitary(U::AbstractMatrix, locs::NTuple{N,Int}) where {N}
    if all(each > 0 for each in tuple_diff(locs))
        return U
    else
        return reorder(U, TupleTools.sortperm(locs))
    end
end

using LinearAlgebra: Transpose
Base.convert(::Type{Transpose{T,Matrix{T}}}, arr::AbstractMatrix{T}) where {T} =
    transpose(Matrix(transpose(arr)))
Base.convert(t::Type{Transpose{T,Matrix{T}}}, arr::Transpose{T}) where {T} =
    invoke(convert, Tuple{Type{Transpose{T,Matrix{T}}},Transpose}, t, arr)
