using Test: @test

using Public

using GSEA

########################################

const P1, P2 = (joinpath(GSEA.P1, st) for st in ("set.json", "data.tsv"))

for pa in readdir(GSEA.P2; join = true)

    if basename(pa) != ".keep"

        rm(pa; recursive = true)

    end

end

const P3, P4, P5, P6, P7 = (
    mkpath(joinpath(GSEA.P2, st)) for st in (
        "data_rank",
        "user_rank",
        "metric_rank.sample",
        "metric_rank.set",
        "user_rank.metric",
    )
)

function read(pa, st = "enrichment")

    Public.make_part(Public.read_table(joinpath(pa, "$st.tsv")))

end

########################################

GSEA.data_rank(P3, P2, P1; minimum = 15, maximum = 500)

const _, S1_, _, N1 = read(P3)

const BO_ = map(nu_ -> any(isfinite, nu_), eachrow(N1))

@test S1_[BO_] == [
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_HYPOXIA",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_KRAS_SIGNALING_DN",
]

@test findmax(N1[BO_, :]) === (0.756249206577638, CartesianIndex(2, 2))

########################################

GSEA.user_rank(P4, joinpath(GSEA.P1, "metric.tsv"), P1)

const _, S2_, _, N2 = read(P4)

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

    @test S2_[nd] === r1

    @test isapprox(N2[nd, :], r2)

end

########################################

const P8 = joinpath(GSEA.P1, "target.tsv")

for (pa, permutation) in ((P5, "sample"), (P6, "set"))

    GSEA.metric_rank(pa, P8, P2, P1; permutation)

end

const (S1, _, S3_, N3), (S2, _, S4_, N4) =
    (read(pa, "metric") for pa in (P5, P6))

@test S1 === S2 === "Gene"

@test S3_[] === S4_[] === "signal-to-noise-ratio"

@test N3 == N4

@test sort!(vec(N4))[[1, 1000]] == [-1.8372355409610066, 1.7411005104346835]

const (_, S5_, _, N5), (_, S6_, _, N6) = (read(pa) for pa in (P5, P6))

@test S5_ == S6_

@test N5[:, 1] == N6[:, 1]

########################################

GSEA.user_rank(P7, joinpath(P6, "metric.tsv"), P1)

const _, S7_, _, N7 = read(P7)

@test S6_ == S7_

@test N6 == N7
