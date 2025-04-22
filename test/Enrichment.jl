using StatsBase: Weights, sample

using Test: @test

using Nucleus

using GSEA

include("_.jl")

# ---- #

# 4.166 ns (0 allocations: 0 bytes)
# 49.291 μs (0 allocations: 0 bytes)
# 3.916 ns (0 allocations: 0 bytes)
# 55.958 μs (0 allocations: 0 bytes)

const A1 = AL_[1]

const A3 = AL_[3]

const NU_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const BO_ = [true, false, true, false, true, true, false, false, true]

const R1_ = randn(100000)

const R2_ = sample([false, true], Weights([0.9, 0.1]), lastindex(R1_))

# TODO: Split

for (al, nu_, bo_, r1, r2) in (
    (A1, NU_, BO_, -0.25, 0.15625),
    (A1, R1_, R2_, nothing, nothing),
    (A3, NU_, BO_, 0.15625, 0.09615384615384615),
    (A3, R1_, R2_, nothing, nothing),
)

    d1, d2 = GSEA.Enrichment.make_delta($al, nu_, bo_)

    @btime GSEA.Enrichment.make_delta($al, $nu_, $bo_)

    @test isnothing(r1) || d1 === r1

    @test isnothing(r2) || d2 === r2

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
# TODO: Split

# 10.553 ns (0 allocations: 0 bytes)
# 10.343 ns (0 allocations: 0 bytes)
# 84.371 ns (0 allocations: 0 bytes)
# 85.883 ns (0 allocations: 0 bytes)
# 154.379 ns (0 allocations: 0 bytes)
# 22.166 μs (0 allocations: 0 bytes)
# 21.541 μs (0 allocations: 0 bytes)
# 127.125 μs (0 allocations: 0 bytes)
# 135.000 μs (0 allocations: 0 bytes)
# 224.208 μs (0 allocations: 0 bytes)

C1_, IN_ = GSEA.Sort.make(C1_, IN_)

GE_, EX_ = GSEA.Sort.make(GE_, EX_)

for (nu_, bo_, re_) in (
        (IN_, Nucleus.Collection.is_in(C1_, C2_), (-0.5, 0, 0, 0, 0)),
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
