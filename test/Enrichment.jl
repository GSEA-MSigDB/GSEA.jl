using StatsBase: Weights, sample

using Test: @test

include("_.jl")

# ---- #

const A1 = AL_[1]

const A3 = AL_[3]

const B1_ = Nucleus.Collection.is_in(C1_, C2_)

const E1 = 1 / 12

const E2 = 1 / 42

const R1_ = randn(100000)

# TODO: Simplify.
const R2_ = sample([false, true], Weights([0.9, 0.1]), lastindex(R1_))

# ---- #

# 5.708 ns (0 allocations: 0 bytes)
# 49.250 μs (0 allocations: 0 bytes)

for (nu_, bo_, r0, r1) in ((IN_, B1_, -1 / 11, E1), (R1_, R2_, nothing, nothing))

    d0, d1 = GSEA.Enrichment.make_delta(A1, nu_, bo_)

    #@btime GSEA.Enrichment.make_delta(A1, $nu_, $bo_)

    @test isnothing(r0) || d0 === r0

    @test isnothing(r1) || d1 === r1

end

# ---- #

# 5.250 ns (0 allocations: 0 bytes)
# 55.958 μs (0 allocations: 0 bytes)

for (nu_, bo_, r1, r2) in ((IN_, B1_, E1, E2), (R1_, R2_, nothing, nothing))

    d1, d2 = GSEA.Enrichment.make_delta(A3, nu_, bo_)

    #@btime GSEA.Enrichment.make_delta(A3, $nu_, $bo_)

    @test isnothing(r1) || d1 === r1

    @test isnothing(r2) || d2 === r2

end

# ---- #

for (d1, d2, re) in ((E1, E2, 1 / 30),)

    @test GSEA.Enrichment.make_delta(d1, d2) === re

end

# ---- #

const P1 = eps()

const P2 = P1 * 2

for (nu, re) in ((-1, P1), (0, P1), (P2, P2))

    @test GSEA.Enrichment.make_eps(nu) === re

end

# ---- #

C1_, IN_ = GSEA.Sort.make(C1_, IN_)

G1_, EX_ = GSEA.Sort.make(G1_, EX_)

const B2_ = Nucleus.Collection.is_in(G1_, G2_)

# ---- #

function test(al, nu_, bo_, re)

    @test isapprox(GSEA.Enrichment.make!(al, nu_, bo_), re; atol = 1e-15)

    #@btime GSEA.Enrichment.make!($al, $nu_, $bo_)

end

# ---- #

# 10.594 ns (0 allocations: 0 bytes)
# 22.125 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((IN_, B1_, -0.5), (EX_, B2_, 0.7651927829281453))

    test(A1, nu_, bo_, re)

end

# ---- #

# 10.375 ns (0 allocations: 0 bytes)
# 21.541 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((IN_, B1_, 0), (EX_, B2_, 0.41482514169516305))

    test(AL_[2], nu_, bo_, re)

end

# ---- #

# 84.501 ns (0 allocations: 0 bytes)
# 127.166 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((IN_, B1_, 0), (EX_, B2_, 1.2297916337424049))

    test(A3, nu_, bo_, re)

end

# ---- #

# 85.927 ns (0 allocations: 0 bytes)
# 135.250 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((IN_, B1_, 0), (EX_, B2_, 1.1161382190540838))

    test(AL_[4], nu_, bo_, re)

end

# ---- #

# 154.637 ns (0 allocations: 0 bytes)
# 224.834 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((IN_, B1_, 0), (EX_, B2_, 1.1181841586127337))

    test(AL_[5], nu_, bo_, re)

end
