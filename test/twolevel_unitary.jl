using YaoArrayRegister: swapcols!, uncols!, mulcol!, unrows!
using Test, LuxurySparse, SparseArrays, LinearAlgebra, StaticArrays

@testset "mulcol!, swapcols!" begin
    a = [1,2,3,4.0]
    mulcol!(a, 2, 0.3)
    @test a ≈ [1, 0.6, 3, 4.0]
    swapcols!(a, 2, 3)
    @test a ≈ [1, 3, 0.6, 4.0]
    swapcols!(a, 2, 3, 0.1, 0.2)
    @test a ≈ [1, 0.12, 0.3, 4.0]

    a = [1 2 3 4.0; 5 6 7 8]
    mulcol!(a, 2, 0.3)
    @test a ≈ [1 0.6 3 4.0; 5 1.8 7 8]
    swapcols!(a, 2, 3)
    @test a ≈ [1 3 0.6 4.0; 5 7 1.8 8]
    swapcols!(a, 2, 3, 0.1, 0.2)
    @test a ≈ [1 0.12 0.3 4.0; 5 0.36 0.7 8]
end

@testset "unrows!" begin
    pm = PermMatrix([2,1], [0.1, 0.2])
    id = IMatrix{2}()
    ds = [1 2; 3.0 4im]
    sp = sparse([1 2; 0 0.1im])
    dg = Diagonal([1 2.0im])
    indices = SVector((2,3))
    for a in [randn(ComplexF64, 4), randn(ComplexF64, 4,4)]
        for m in [pm, id, ds, dg]
            for m_ in [m, staticize(m)]
                a1 = unrows!(copy(a), indices, m_)
                ta = ndims(a) == 2 ? copy(transpose(a)) : a
                a2 = uncols!(ta, indices, m_)
                a2 = ndims(a2) == 2 ? transpose(a2) : a2
                @test a1 ≈ a2
                if m_ !== m
                    @show typeof(m_)
                    @test (@allocated unrows!(a1, indices, m_)) == 0
                    @test (@allocated uncols!(a2, indices, m_)) == 0
                end
            end
        end
    end
end
