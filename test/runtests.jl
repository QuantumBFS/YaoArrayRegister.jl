using YaoArrayRegister
using Test
using Aqua
Aqua.test_all(YaoArrayRegister)

@testset "test utils" begin
    include("utils.jl")
end

@testset "test ArrayReg" begin
    include("register.jl")
    include("operations.jl")
    include("focus.jl")
end

@testset "test instructions" begin
    include("instruct.jl")
end

@testset "test measure" begin
    include("measure.jl")
end
