using PkgBenchmark, BenchmarkTools
using YaoArrayRegister, BitBasis, Random, YaoBase

bench(n, U, loc) = @benchmarkable instruct!(st, $U, $loc) setup=(st=statevec(rand_state($n)))
bench(n, U, loc, control_locs, control_bits) = @benchmarkable instruct!(st, $U, $loc, $control_locs, $control_bits) setup=(st=statevec(rand_state($n)))

const SUITE = BenchmarkGroup()

# Specialized Gate Instruction
SUITE["Single Qubit"] = BenchmarkGroup()
## single qubit benchmark
for U in YaoArrayRegister.SPECIALIZATION_LIST
    println(U)
    for n in 1:4:25
        SUITE["Single Qubit"][string(U)] = @benchmarkable bench($n, Val($(QuoteNode(U))), 1)
    end
end

SUITE["Multi Qubit"] = BenchmarkGroup()
const location_sparsity = 0.4
## multi qubit benchmark
for U in YaoArrayRegister.SPECIALIZATION_LIST
    for n in 1:4:25
        perms = randperm(n)[1:round(Int, location_sparsity * n)]
        SUITE["Multi Qubit"][string(U)] = @benchmarkable bench($n, Val($(QuoteNode(U))), $perms)
    end
end

# General Instructions
SUITE["General Matrix"] = BenchmarkGroup()
## General Matrix Instruction

for n in 1:10, T in [ComplexF32, ComplexF64]
    for U in [
        rand_unitary(T, 1<<n), # dense matrices
        sprand_unitary(T, 1<<n),
        sprand_hermitian(T, 1<<n),
        ]
    end
end

## Special Matrices Instruction
SUITE["Special Matrices"] = BenchmarkGroup()
