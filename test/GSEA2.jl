using Test: @test

using Nucleus

using GSEA

include("_.jl")

# ---- #

const T1 = joinpath(DA, "target.tsv")

const T2 = joinpath(DA, "data.tsv")

const T3 = joinpath(DA, "metric.tsv")

const JS = joinpath(DA, "set.json")

const O1, O2, O3, O4, O5 =
    (mkpath(joinpath(TE, ba)) for ba in ("da", "us", "me", "us2", "me2"))

const T4 = joinpath(O3, "metric.tsv")

function rea(di)

    Nucleus.Table.rea(joinpath(di, "result.tsv"))

end

# ---- #
# TODO

GSEA.data_rank(O1, T2, JS; minimum = 15, maximum = 500)

const A1 = rea(O1)

@test size(A1) === (8, 10)

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
# TODO

GSEA.user_rank(O2, T3, JS)

const A2 = rea(O2)

@test size(A2) === (50, 5)

for (id, r1, r2) in (
    (
        45,
        "HALLMARK_PANCREAS_BETA_CELLS",
        [-0.3526604388911228, -1.424806498309897, 0.01718494271685761, 0.13747954173486088],
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
        [0.8665786760826208, 3.27574664008412, 0.00026469031233456857, 0.00277924827951297],
    ),
)

    @test A2[id, 1] === r1

    @test is_egal(Vector(A2[id, 2:end]), r2)

end

# ---- #
# TODO

GSEA.metric_rank(O3, T1, T2, JS)

const A4 = Nucleus.Table.rea(T4)

@test size(A4) === (1000, 2)

@test is_egal(names(A4), ["Feature", "signal-to-noise-ratio"])

@test is_egal(sort(A4, 2)[[1, end], 2], [-1.8372355409610066, 1.7411005104346835])

const A3 = rea(O3)

@test size(A3) === (50, 5)

# ---- #

GSEA.user_rank(O4, T4, JS)

const A5 = rea(O4)

@test is_egal(A5[!, 1:2], A3[!, 1:2])

# ---- #

GSEA.metric_rank(O5, T1, T2, JS; permutation = "set")

@test is_egal(rea(O5), A5)
