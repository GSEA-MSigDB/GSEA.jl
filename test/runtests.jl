using GSEA

# ------------------------------------ #

using BenchmarkTools: @btime

using Test: @test

using Public

using GSEA

########################################

const AL_ = GSEA.S0(), GSEA.S0a(), GSEA.D2(), GSEA.D2f(), GSEA.D0f2f()

########################################

const S1_, N1_ = GSEA.make_sort(
    eachcol(Public.read_table(joinpath(GSEA.P1, "myc.tsv"))[!, 1:2])...,
)

const B1_ = map(
    in(
        convert(
            Dict{String, Vector{String}},
            Public.read_pair(joinpath(GSEA.P1, "c2.all.v7.1.symbols.json")),
        )["COLLER_MYC_TARGETS_UP"],
    ),
    S1_,
)

########################################

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

    @test GSEA.number_enrichment!(al, N1_, B1_) === re

    @btime GSEA.number_enrichment!($al, N1_, B1_)

end

########################################

const P1 = joinpath(GSEA.P1, "h.all.v7.1.symbols.json")

########################################

const S2_, ST__ = GSEA.read_pair(P1)

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

for st_ in (['A', 'I', '_'], ['D', 'E', 'F', '_']), al in AL_

    GSEA.write_enrichment(
        "",
        al,
        'A':'I',
        -4:4.0,
        st_,
        Dict(Public.pair_title("$al")),
    )

end

########################################

const P2 = joinpath(GSEA.P1, "data.tsv")

for pa in readdir(GSEA.P2; join = true)

    if basename(pa) != ".keep"

        rm(pa; recursive = true)

    end

end

const P3, P4, P5, P6, P7 = (
    mkpath(joinpath(GSEA.P2, st)) for st in (
        "data_rank",
        "user_rank_1",
        "metric_rank_sample",
        "metric_rank_set",
        "user_rank_2",
    )
)

########################################

GSEA.data_rank(P3, P2, P1; minimum = 15, maximum = 500)

const _, S3_, _, N1 =
    Public.make_part(Public.read_table(joinpath(P3, "result.tsv")))

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

const _, S4_, _, N2 =
    Public.make_part(Public.read_table(joinpath(P4, "result.tsv")))

for (nd, r1, r2) in (
    (
        1,
        "HALLMARK_ADIPOGENESIS",
        [
            0.3285678436753619,
            1.5493934448392406,
            0.005675675675675676,
            0.02648648648648,
        ],
        649,
    ),
    (
        50,
        "HALLMARK_XENOBIOTIC_METABOLISM",
        [
            0.2515741729236317,
            1.173484194980481,
            0.17108108108108108,
            0.2993918918918919,
        ],
    ),
)

    @test S4_[nd] === r1

    @test isapprox(N2[nd, :], r2)

end

########################################

GSEA.metric_rank(
    P5,
    joinpath(GSEA.P1, "target.tsv"),
    P2,
    P1;
    permutation = "sample",
)

GSEA.metric_rank(
    P6,
    joinpath(GSEA.P1, "target.tsv"),
    P2,
    P1;
    permutation = "set",
)

const (S1, _, S5_, N3), (S2, _, S6_, N4) = (
    Public.make_part(Public.read_table(joinpath(pa, "metric.tsv"))) for
    pa in (P5, P6)
)

@test S1 === S2 === "Gene"

@test S5_[] === S6_[] === "signal-to-noise-ratio"

@test N3 == N4

@test sort!(vec(N4))[[1, 1000]] == [-1.8372355409610066, 1.7411005104346835]

const (_, S7_, _, N5), (_, S8_, _, N6) = (
    Public.make_part(Public.read_table(joinpath(pa, "result.tsv"))) for
    pa in (P5, P6)
)

@test S7_ == S8_

@test N5[:, 1] == N6[:, 1]

########################################

GSEA.user_rank(P7, joinpath(P6, "metric.tsv"), P1)

const _, S9_, _, N7 =
    Public.make_part(Public.read_table(joinpath(P7, "result.tsv")))

@test S8_ == S9_

@test N6 == N7
