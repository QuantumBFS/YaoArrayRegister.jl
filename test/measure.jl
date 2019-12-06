using Test, YaoArrayRegister, YaoBase

@testset "select" begin
    reg = product_state(4, 6; nbatch = 2)
    # println(focus!(reg, [1,3]))
    r1 = select!(focus!(copy(reg), [2, 3]), 0b11) |> relax!(to_nactive = 2)
    r2 = select(focus!(copy(reg), [2, 3]), 0b11) |> relax!(to_nactive = 2)
    r3 = copy(reg) |> focus!(2, 3) |> select!(0b11) |> relax!(to_nactive = 2)

    @test r1' * r1 ≈ ones(2)
    @test r1 ≈ r2
    @test r3 ≈ r2
end

@testset "measure and resetto/remove" begin
    reg = rand_state(4)
    res = measure_resetto!(reg, (4,))
    @test isnormalized(reg)
    result = measure(reg; nshots = 10)
    @test all(result .< 8)

    reg = rand_state(6) |> focus!(1, 4, 3)
    reg0 = copy(reg)
    res = measure_remove!(reg)
    select(reg0, res)
    @test select(reg0, res) |> normalize! ≈ reg

    r = rand_state(10)
    r1 = copy(r) |> focus!(1, 4, 3)
    res = measure_remove!(r, (1, 4, 3))
    r2 = select(r1, res)
    r2 = relax!(r2, (); to_nactive = nqubits(r2))
    @test normalize!(r2) ≈ r

    reg = rand_state(6, nbatch = 5) |> focus!((1:5)...)
    measure_resetto!(reg, 1)
    @test nactive(reg) == 5
end
