using BenchmarkTools: @btime

using Test: @test

using Public

using GSEA

########################################

const AL_ = GSEA.KS0(), GSEA.A0(), GSEA.DA2(), GSEA.DA2W(), GSEA.DA2W0W()

########################################

const N1_ = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6.0]

const B1_ = [
    true,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    true,
]

const N1 = 1 / 12

const N2 = 1 / 42

for (nd, nu_, bo_, r1, r2) in
    ((1, N1_, B1_, -1 / 11, N1), (3, N1_, B1_, N1, N2))

    GSEA.number_delta(AL_[nd], nu_, bo_) === (r1, r2)

end

@test GSEA.number_delta(N1, N2) === 1 / 30

########################################

const EP = eps()

for (nu, re) in ((-1, EP), (0, EP), fill(EP * 2, 2))

    @test GSEA.number_eps(nu) === re

end

########################################

const S1_, N2_ = GSEA.make_sort(
    eachcol(Public.read_table(joinpath(GSEA.P1, "myc.tsv"))[!, 1:2])...,
)

########################################

const B2_ = map(
    in(
        GSEA.read_gmt(joinpath(GSEA.P1, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"],
    ),
    S1_,
)

# 22.208 μs (0 allocations: 0 bytes)
# 21.583 μs (0 allocations: 0 bytes)
# 123.292 μs (0 allocations: 0 bytes)
# 134.542 μs (0 allocations: 0 bytes)
# 218.792 μs (0 allocations: 0 bytes)
for (nd, nu_, bo_, re) in (
    (1, N2_, B2_, 0.7651927829281453),
    (2, N2_, B2_, 0.41482514169516305),
    (3, N2_, B2_, 1.2297916337424049),
    (4, N2_, B2_, 1.1161382190540838),
    (5, N2_, B2_, 1.1181841586127337),
)

    al = AL_[nd]

    @test GSEA.number_enrichment!(al, nu_, bo_) === re

    @btime GSEA.number_enrichment!($al, $nu_, $bo_)

end

########################################

const DI = GSEA.read_gmt(joinpath(GSEA.P1, "h.all.v7.1.symbols.gmt"))

const S2_ = collect(keys(DI))

const ST__ = collect(values(DI))

# 1.622 ms (32 allocations: 1.88 MiB)
# 1.619 ms (32 allocations: 1.88 MiB)
# 7.284 ms (32 allocations: 1.88 MiB)
# 7.815 ms (32 allocations: 1.88 MiB)
# 12.704 ms (32 allocations: 1.88 MiB)
for al in AL_

    @test S2_[sortperm(GSEA.number_enrichment(al, S1_, N2_, ST__))][49:50] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

    @btime GSEA.number_enrichment($al, S1_, N2_, ST__)

end
