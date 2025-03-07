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
# 86.625 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((N1_, B1_, (-0.25, 0.15625)), (N2_, B2_, nothing))

    al = AL_[1]

    @test isnothing(re) || GSEA.Enrichment.make_delta(al, nu_, bo_) === re

    #@btime GSEA.Enrichment.make_delta($al, $nu_, $bo_)

end

# ---- #

# 6.417 ns (0 allocations: 0 bytes)
# 92.833 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((N1_, B1_, (0.15625, 0.09615384615384615)), (N2_, B2_, nothing))

    al = AL_[3]

    @test isnothing(re) || GSEA.Enrichment.make_delta(al, nu_, bo_) === re

    #@btime GSEA.Enrichment.make_delta($al, $nu_, $bo_)

end

# ---- #

for (d1, d2, re) in ((1 / 3, 0.5, -1.0),)

    @test GSEA.Enrichment.make_delta(d1, d2) === re

end

# ---- #

const EP = eps()

for (po, re) in ((-EP, EP), (0, EP))

    @test GSEA.Enrichment.make_eps(po) === re

end

# ---- #

# 15.708 ns (0 allocations: 0 bytes)
# 15.698 ns (0 allocations: 0 bytes)
# 154.755 ns (0 allocations: 0 bytes)
# 155.930 ns (0 allocations: 0 bytes)
# 284.304 ns (0 allocations: 0 bytes)
# 45.416 μs (0 allocations: 0 bytes)
# 37.541 μs (0 allocations: 0 bytes)
# 224.542 μs (0 allocations: 0 bytes)
# 243.916 μs (0 allocations: 0 bytes)
# 414.875 μs (0 allocations: 0 bytes)

C1_, LI_ = GSEA.Sort.make(C1_, LI_)

GE_, EX_ = GSEA.Sort.make(GE_, EX_)

for (nu_, bo_, re_) in (
        (LI_, Nucleus.Collection.is_in(C1_, C2_), (-0.5, 0.0, 0.0, 0.0, 0.0)),
        (
            EX_,
            Nucleus.Collection.is_in(GE_, D2["COLLER_MYC_TARGETS_UP"]),
            (
                0.7651927829281453,
                0.41482514169516305,
                1.2297916337424049,
                1.1161382190540838,
                1.1181841586127337,
            ),
        ),
    ),
    (al, re) in zip(AL_, re_)

    @test isapprox(GSEA.Enrichment.make!(al, nu_, bo_), re; atol = 1e-15)

    #@btime GSEA.Enrichment.make!($al, $nu_, $bo_)

end
