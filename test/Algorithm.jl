using Test: @test

using GSEA

using StatsBase: Weights, sample

using Nucleus

include("_.jl")

# ---- #

const N1_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const B1_ = [true, false, true, false, true, true, false, false, true]

const N2_ = randn(100000)

const B2_ = sample([false, true], Weights([0.9, 0.1]), lastindex(N2_))

# ---- #

# 7.958 ns (0 allocations: 0 bytes)
# 86.500 μs (0 allocations: 0 bytes)

const NO = -0.25

for (nu_, bo_, re) in ((N1_, B1_, (NO, 0.15625)), (N2_, B2_, nothing))

    al = GSEA.Algorithm.KS()

    @test isnothing(re) || GSEA.Algorithm.make_normalizer(al, nu_, bo_) === re

    #@btime GSEA.Algorithm.make_normalizer($al, $nu_, $bo_)

end

# ---- #

# 6.417 ns (0 allocations: 0 bytes)
# 92.833 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((N1_, B1_, (0.15625, 0.09615384615384615)), (N2_, B2_, nothing))

    al = GSEA.Algorithm.KLioM()

    @test isnothing(re) || GSEA.Algorithm.make_normalizer(al, nu_, bo_) === re

    #@btime GSEA.Algorithm.make_normalizer($al, $nu_, $bo_)

end

# ---- #

for (o1, o2, re) in ((1 / 3, 0.5, -1.0),)

    @test GSEA.Algorithm.make_normalizer(o1, o2) === re

end

# ---- #

const EP = eps()

for (nu, re) in ((-EP, EP), (0, EP))

    @test GSEA.Algorithm.make_eps(nu) === re

end

# ---- #

# 15.739 ns (0 allocations: 0 bytes)
# 15.698 ns (0 allocations: 0 bytes)
# 281.250 ns (0 allocations: 0 bytes)
# 281.250 ns (0 allocations: 0 bytes)
# 156.030 ns (0 allocations: 0 bytes)
# 154.908 ns (0 allocations: 0 bytes)
# 45.334 μs (0 allocations: 0 bytes)
# 37.500 μs (0 allocations: 0 bytes)
# 410.334 μs (0 allocations: 0 bytes)
# 410.416 μs (0 allocations: 0 bytes)
# 245.375 μs (0 allocations: 0 bytes)
# 224.541 μs (0 allocations: 0 bytes)

C1_, LI_ = GSEA.update(C1_, LI_)

GE_, EX_ = GSEA.update(GE_, EX_)

for (nu_, bo_, re_) in (
        (LI_, Nucleus.Collection.is_in(C1_, C2_), (-0.5, 0.0, 0.0, 0.0, 0.0, 0.0)),
        (
            EX_,
            Nucleus.Collection.is_in(GE_, D2["COLLER_MYC_TARGETS_UP"]),
            (
                0.7651927829281453,
                0.41482514169516305,
                1.1181841586127337,
                1.1140922794954267,
                1.1161382190540838,
                1.2297916337424049,
            ),
        ),
    ),
    (al, re) in zip(AL_, re_)

    @test isapprox(GSEA.Algorithm.make!(al, nu_, bo_, nothing), re; atol = 1e-15)

    #@btime GSEA.Algorithm.make!($al, $nu_, $bo_, nothing)

end
