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

########################################

const B1_ = map(
    in(
        convert(
            Dict{String, Vector{String}},
            Public.read_pair(joinpath(GSEA.P1, "c2.all.v7.1.symbols.json")),
        )["COLLER_MYC_TARGETS_UP"],
    ),
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

    @test GSEA.number_enrichment!(al, N1_, B1_) === re

    @btime GSEA.number_enrichment!($al, N1_, B1_)

end

########################################

const P1 = joinpath(GSEA.P1, "h.all.v7.1.symbols.json")

const DI::Dict{String, Vector{String}} = Public.read_pair(P1)

const S2_ = collect(keys(DI))

const ST__ = collect(values(DI))

# 1.618 ms (32 allocations: 1.88 MiB)
# 1.606 ms (32 allocations: 1.88 MiB)
# 7.232 ms (32 allocations: 1.88 MiB)
# 7.629 ms (32 allocations: 1.88 MiB)
# 12.515 ms (32 allocations: 1.88 MiB)
for al in AL_

    @test S2_[sortperm(GSEA.number_enrichment(al, S1_, N1_, ST__))][49:50] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

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
        Dict("title" => Dict("text" => "$al")),
    )

end

########################################

const P2 = cp(pkgdir(GSEA, "ou"), joinpath(tempdir(), "GSEA"); force = true)

const P3 = joinpath(GSEA.P1, "data.tsv")

########################################

const P4 = mkpath(joinpath(P2, "data_rank"))

GSEA.data_rank(P4, P3, P1; minimum = 15, maximum = 500)

const _, S1_, _, N1 =
    Public.make_part(Public.read_table(joinpath(P4, "result.tsv")))

const B2_ = map(nu_ -> !any(isnan, nu_), eachrow(N1))

@test S1_[B2_] == [
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
]

@test findmax(N1[B2_, :]) === (0.756249206577638, CartesianIndex(2, 2))

########################################

const P5 = mkpath(joinpath(P2, "user_rank_1"))

GSEA.user_rank(P5, joinpath(GSEA.P1, "metric.tsv"), P1)

const _, S2_, _, N2 =
    Public.make_part(Public.read_table(joinpath(P5, "result.tsv")))

for (nd, r1, r2) in (
    (
        45,
        "HALLMARK_PANCREAS_BETA_CELLS",
        [
            -0.3526604388911228,
            -1.424806498309897,
            0.01718494271685761,
            0.13747954173486088,
        ],
    ),
    (
        33,
        "HALLMARK_PROTEIN_SECRETION",
        [
            -0.27209592258934645,
            -1.324926701473525,
            0.03764320785597381,
            0.15057283142389524,
        ],
    ),
    (
        36,
        "HALLMARK_MYC_TARGETS_V1",
        [
            0.6033555857451218,
            2.719800888801976,
            0.00026469031233456857,
            0.00277924827951297,
        ],
    ),
    (
        10,
        "HALLMARK_MYC_TARGETS_V2",
        [
            0.8665786760826208,
            3.27574664008412,
            0.00026469031233456857,
            0.00277924827951297,
        ],
    ),
)

    @test S2_[nd] === r1

    @test isapprox(N2[nd, :], r2)

end

########################################

const S8_ = "sample", "set"

const PA_ = map(st -> mkpath(joinpath(P2, "metric_rank_$st")), S8_)

for nd in 1:2

    GSEA.metric_rank(
        PA_[nd],
        joinpath(GSEA.P1, "target.tsv"),
        P3,
        P1;
        permutation = S8_[nd],
    )

end

const S1, _, S3_, N3 =
    Public.make_part(Public.read_table(joinpath(PA_[1], "metric.tsv")))

const S2, _, S4_, N4 =
    Public.make_part(Public.read_table(joinpath(PA_[2], "metric.tsv")))

@test S1 === S2 === "Gene"

@test S3_[] === S4_[] === "signal-to-noise-ratio"

@test N3 == N4

@test sort!(vec(N3))[[1, 1000]] == [-1.8372355409610066, 1.7411005104346835]

const _, S5_, _, N5 =
    Public.make_part(Public.read_table(joinpath(PA_[1], "result.tsv")))

const _, S6_, _, N6 =
    Public.make_part(Public.read_table(joinpath(PA_[2], "result.tsv")))

@test S5_ == S6_

@test N5[:, 1] == N6[:, 1]

########################################

const P6 = mkpath(joinpath(P2, "user_rank_2"))

GSEA.user_rank(P6, joinpath(PA_[2], "metric.tsv"), P1)

const _, S7_, _, N7 =
    Public.make_part(Public.read_table(joinpath(P6, "result.tsv")))

@test S6_ == S7_

@test N6 == N7
