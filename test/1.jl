using BenchmarkTools: @btime

using Test: @test

using Public

using GSEA

########################################

const AL_ = GSEA.S0(), GSEA.S0a(), GSEA.D2(), GSEA.D2f(), GSEA.D0f2f()

########################################

const S1_, N1_ = GSEA.make_sort(
    eachcol(Public.read_table(joinpath(GSEA.P1, "metric.tsv"))[!, 1:2])...,
)

########################################

const BO_ = map(
    in((
        "AHCY",
        "AK4",
        "ASS1",
        "C1QBP",
        "CCND2",
        "CEBPZ",
        "CKS2",
        "EIF5A",
        "FABP5",
        "FBL",
        "G0S2",
        "GPI",
        "GRPEL1",
        "HDGF",
        "HSPD1",
        "IARS1",
        "NAMPT",
        "NCL",
        "ODC1",
        "POLR2H",
        "PPIF",
        "SLC16A1",
        "TFRC",
        "TRAP1",
    )),
    S1_,
)

# 22.167 μs (0 allocations: 0 bytes)
# 21.542 μs (0 allocations: 0 bytes)
# 123.458 μs (0 allocations: 0 bytes)
# 133.334 μs (0 allocations: 0 bytes)
# 217.000 μs (0 allocations: 0 bytes)
for (nd, re) in (
    (1, 0.7651927829281453),
    (2, 0.41482514169516305),
    (3, 1.2297916337424049),
    (4, 1.1161382190540838),
    (5, 1.1181841586127337),
)

    al = AL_[nd]

    @test GSEA.number_enrichment!(al, N1_, BO_) === re

    @btime GSEA.number_enrichment!($al, N1_, BO_)

end

########################################

const S2_, ST__ = GSEA.read_pair(joinpath(GSEA.P1, "set.json"))

# 1.618 ms (32 allocations: 1.88 MiB)
# 1.606 ms (32 allocations: 1.88 MiB)
# 7.232 ms (32 allocations: 1.88 MiB)
# 7.629 ms (32 allocations: 1.88 MiB)
# 12.515 ms (32 allocations: 1.88 MiB)
for al in AL_

    @test S2_[partialsortperm(
        GSEA.number_enrichment(al, S1_, N1_, ST__),
        49:50,
    )] == ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

    @btime GSEA.number_enrichment($al, S1_, N1_, ST__)

end

########################################

const CH_ = 'A':'I'

const IN_ = -4:4

for ch_ in (['_', 'A', 'I'], ['D', 'E', 'F', '_']), al in AL_

    GSEA.write_enrichment("", al, CH_, IN_, ch_, Dict(Public.pair_title(al)))

end
