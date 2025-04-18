using Test: @test

using GSEA

using Nucleus

include("_.jl")

# ---- #

const JS = joinpath(DA, "set.json")

const T1 = joinpath(DA, "data.tsv")

# ---- #

const O1 = mkpath(joinpath(TE, "data_rank"))

GSEA.data_rank(O1, T1, JS; minimum = 15, maximum = 500)

const A1 = Nucleus.Table.rea(joinpath(O1, "result.tsv"))

@test is_egal(
    A1[!, 1],
    [
        "HALLMARK_ESTROGEN_RESPONSE_LATE",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_ESTROGEN_RESPONSE_EARLY",
        "HALLMARK_KRAS_SIGNALING_DN",
        "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_APICAL_JUNCTION",
        "HALLMARK_HYPOXIA",
        "HALLMARK_GLYCOLYSIS",
    ],
)

@test findmax(Matrix(A1[!, 2:end])) === (0.756249206577638, CartesianIndex(2, 2))

# ---- #

const KE_ = (minimum = 15, maximum = 500, more_plots = "HALLMARK_UV_RESPONSE_DN")

# ---- #

const O2 = mkpath(joinpath(TE, "user_rank"))

GSEA.user_rank(
    O2,
    joinpath(DA, "metric.tsv"),
    JS;
    more_plots = "HALLMARK_MYC_TARGETS_V1;HALLMARK_UV_RESPONSE_DN;HALLMARK_UV_RESPONSE_UP;ALIEN",
)

const A2 = Nucleus.Table.rea(joinpath(O2, "result.tsv"))

@test size(A2, 1) === 50

for (id, r1, r2, r3, r4) in (
    (45, "HALLMARK_PANCREAS_BETA_CELLS", -0.35266, -1.36616, 0.0200837),
    (33, "HALLMARK_PROTEIN_SECRETION", -0.272096, -1.25207, 0.0686192),
    (36, "HALLMARK_MYC_TARGETS_V1", 0.603356, 2.73998, 0.000262812),
    (10, "HALLMARK_MYC_TARGETS_V2", 0.866579, 3.36557, 0.000262812),
)

    @test A2[id, 1] === r1

    @test isapprox(A2[id, 2], r2; atol = 1e-6)

    # TODO

    #@test isapprox(A2[id, 3], r3; atol = 1e-5)

    #@test isapprox(A2[id, 4], r4; atol = 1e-7)

end

# ---- #

const O3 = mkpath(joinpath(TE, "metric_rank"))

const T2 = joinpath(DA, "target.tsv")

const T3 = joinpath(O3, "metric.tsv")

# ---- #

GSEA.metric_rank(O3, T2, T1, JS; KE_...)

const A3 = Nucleus.Table.rea(T3)

const A4 = Nucleus.Table.rea(joinpath(O3, "result.tsv"))

@test size(A3) === (1000, 2)

@test is_egal(names(A3), ["Feature", "signal-to-noise-ratio"])

@test is_egal(sort(A3, 2)[[1, end], 2], [-1.8372355409610066, 1.7411005104346835])

@test size(A4, 1) === 8

# ---- #

GSEA.user_rank(TE, T3, JS; KE_...)

const A5 = Nucleus.Table.rea(joinpath(TE, "result.tsv"))

@test is_egal(A5[!, 1:2], A4[!, 1:2])

# ---- #

GSEA.metric_rank(TE, T2, T1, JS; permutation = "set", KE_...)

const A6 = Nucleus.Table.rea(joinpath(TE, "result.tsv"))

@test is_egal(A6, A5)
