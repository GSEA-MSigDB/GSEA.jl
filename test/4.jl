using Test: @test

using Public

using GSEA

########################################

const PA = cp(pkgdir(GSEA, "ou"), joinpath(tempdir(), "GSEA"); force = true)

########################################

GSEA.data_rank(
    mkpath(joinpath(PA, "data_rank")),
    joinpath(GSEA.P1, "data.tsv"),
    joinpath(GSEA.P1, "set.json");
    minimum = 15,
    maximum = 500,
)

const _, S1_, _, N1 =
    Public.make_part(Public.read_table(joinpath(PA, "data_rank", "result.tsv")))

const BO_ = map(nu_ -> !any(isnan, nu_), eachrow(N1))

@test S1_[BO_] == [
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
]

@test findmax(N1[BO_, :]) === (0.756249206577638, CartesianIndex(2, 2))

########################################

GSEA.user_rank(
    mkpath(joinpath(PA, "user_rank_1")),
    joinpath(GSEA.P1, "metric.tsv"),
    joinpath(GSEA.P1, "set.json"),
)

const _, S2_, _, N2 = Public.make_part(
    Public.read_table(joinpath(PA, "user_rank_1", "result.tsv")),
)

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

    @test N2[nd, :] == r2

end

########################################

for permutation in ("sample", "set")

    GSEA.metric_rank(
        mkpath(joinpath(PA, "metric_rank_$permutation")),
        joinpath(GSEA.P1, "target.tsv"),
        joinpath(GSEA.P1, "data.tsv"),
        joinpath(GSEA.P1, "set.json");
        permutation,
    )

end

const S1, _, S3_, N3 = Public.make_part(
    Public.read_table(joinpath(PA, "metric_rank_sample", "metric.tsv")),
)

const S2, _, S4_, N4 = Public.make_part(
    Public.read_table(joinpath(PA, "metric_rank_set", "metric.tsv")),
)

@test S1 === S2 === "Gene"

@test S3_[] === S4_[] === "signal-to-noise-ratio"

@test N3 == N4

@test sort!(vec(N3))[[1, 1000]] == [-1.8372355409610066, 1.7411005104346835]

const _, S5_, _, N5 = Public.make_part(
    Public.read_table(joinpath(PA, "metric_rank_sample", "result.tsv")),
)

const _, S6_, _, N6 = Public.make_part(
    Public.read_table(joinpath(PA, "metric_rank_set", "result.tsv")),
)

@test S5_ == S6_

@test N5[:, 1] == N6[:, 1]

########################################

GSEA.user_rank(
    mkpath(joinpath(PA, "user_rank_2")),
    joinpath(PA, "metric_rank_set", "metric.tsv"),
    joinpath(GSEA.P1, "set.json"),
)

const _, S7_, _, N7 = Public.make_part(
    Public.read_table(joinpath(PA, "user_rank_2", "result.tsv")),
)

@test S6_ == S7_

@test N6 == N7
