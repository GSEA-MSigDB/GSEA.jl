using GSEA

# ------------------------------------ #

using BenchmarkTools: @btime

using Test: @test

using Public

########################################

const AL_ = GSEA.S0(), GSEA.S0a(), GSEA.D2(), GSEA.D2w(), GSEA.DD()

########################################

const S1_, N1_ = GSEA.make_score(
    eachcol(Public.read_table(joinpath(GSEA.P1, "metric.tsv"))[!, 1:2])...,
)

########################################

const B1_ = map(
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

# 19.666 μs (0 allocations: 0 bytes)
# 19.875 μs (0 allocations: 0 bytes)
# 132.625 μs (0 allocations: 0 bytes)
# 144.834 μs (0 allocations: 0 bytes)
# 231.334 μs (0 allocations: 0 bytes)
for (nd, re) in (
    (1, 0.7651927829281453),
    (2, 0.41482514169516305),
    (3, 1.2297916337424049),
    (4, 1.1161382190540838),
    (5, 1.1181841586127337),
)

    al = AL_[nd]

    @test GSEA.number_enrichment!(al, N1_, B1_) === re

    @btime GSEA.number_enrichment!($al, N1_, B1_)

end

########################################

for ch_ in (['_', 'A', 'I'], ['D', 'E', 'F', '_']), al in AL_

    GSEA.write_enrichment(
        "",
        al,
        'A':'I',
        -4:4,
        ch_,
        Dict(Public.pair_title(al)),
    )

end

########################################

const S2_, ST__ = GSEA.read_pair(joinpath(GSEA.P1, "set.json"))

# 2.087 ms (31 allocations: 1.90 MiB)
# 2.088 ms (31 allocations: 1.90 MiB)
# 8.311 ms (31 allocations: 1.90 MiB)
# 8.775 ms (31 allocations: 1.90 MiB)
# 13.399 ms (31 allocations: 1.90 MiB)
for al in AL_

    @test S2_[partialsortperm(
        GSEA.number_enrichment(al, S1_, N1_, ST__),
        49:50,
    )] == ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

    @btime GSEA.number_enrichment($al, S1_, N1_, ST__)

end

########################################

const P1 = joinpath(GSEA.P1, "set.json")

const P2 = joinpath(GSEA.P1, "data.tsv")

for st in filter!(!=(".keep"), readdir(GSEA.P2))

    rm(joinpath(GSEA.P2, st); recursive = true)

end

const P3 = mkpath(joinpath(GSEA.P2, "data_rank"))

const P4 = mkpath(joinpath(GSEA.P2, "user_rank"))

const P5 = mkpath(joinpath(GSEA.P2, "metric_rank.sample"))

const P6 = mkpath(joinpath(GSEA.P2, "metric_rank.set"))

const P7 = mkpath(joinpath(GSEA.P2, "user_rank.metric"))

function read(pa, st = "enrichment")

    Public.make_part(Public.read_table(joinpath(pa, "$st.tsv")))

end

########################################

GSEA.data_rank(P3, P2, P1; minimum = 15, maximum = 500)

const _, S3_, _, N1 = read(P3)

const B2_ = map(nu_ -> any(isfinite, nu_), eachrow(N1))

@test S3_[B2_] == [
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_HYPOXIA",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_KRAS_SIGNALING_DN",
]

@test findmax(N1[B2_, :]) === (0.756249206577638, CartesianIndex(2, 2))

########################################

GSEA.user_rank(P4, joinpath(GSEA.P1, "metric.tsv"), P1)

const _, S4_, _, N2 = read(P4)

for (nd, r1, r2) in (
    (
        1,
        "HALLMARK_ADIPOGENESIS",
        [
            0.2567329280012404,
            1.6672452839970535,
            0.005025125628140704,
            0.0076005025125628145,
        ],
        649,
    ),
    (
        50,
        "HALLMARK_XENOBIOTIC_METABOLISM",
        [
            -0.22069178947934942,
            -1.2891137391302299,
            0.08372093023255814,
            0.09894291754756872,
        ],
    ),
)

    @test S4_[nd] === r1

    @test isapprox(N2[nd, :], r2)

end

########################################

const P8 = joinpath(GSEA.P1, "target.tsv")

for (pa, permutation) in ((P5, "sample"), (P6, "set"))

    GSEA.metric_rank(pa, P8, P2, P1; permutation)

end

const S1, _, S5_, N3 = read(P5, "metric")

const S2, _, S6_, N4 = read(P6, "metric")

@test S1 === S2 === "Gene"

@test S5_[] === S6_[] === "signal-to-noise-ratio"

@test N3 == N4

@test sort!(vec(N4))[[1, 1000]] == [-1.8372355409610066, 1.7411005104346835]

const _, S7_, _, N5 = read(P5)

const _, S8_, _, N6 = read(P6)

@test S7_ == S8_

@test N5[:, 1] == N6[:, 1]

########################################

GSEA.user_rank(P7, joinpath(P6, "metric.tsv"), P1)

const _, S9_, _, N7 = read(P7)

@test S8_ == S9_

@test N6 == N7
