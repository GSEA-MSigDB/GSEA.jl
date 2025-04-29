using Test: @test

include("_.jl")

# ---- #

function test(fu, ba, re)

    ts = joinpath(TE, "_.tsv")

    fu(ts, joinpath(DA, ba))

    @test size(Nucleus.Table.rea(ts)) === re

end

# ---- #

for (ba, re) in (
    ("1.cls", (1, 7)),
    ("GSE76137.cls", (1, 7)),
    ("CCLE_mRNA_20Q2_no_haem_phen.cls", (1, 900)),
)

    test(GSEA.cls, ba, re)

end

# ---- #

for (ba, re) in (("1.gct", (13321, 190)),)

    test(GSEA.gct, ba, re)

end

# ---- #

const UM = 50

# ---- #

const J1 = joinpath(TE, "_.json")

for (ba, re) in (("h.all.v7.1.symbols.gmt", UM), ("c2.all.v7.1.symbols.gmt", 5529))

    GSEA.gmt(J1, joinpath(DA, ba))

    @test length(Nucleus.Dictionary.rea(J1)) === re

end

# ---- #

const T1 = joinpath(DA, "target.tsv")

const T2 = joinpath(DA, "data.tsv")

const T3 = joinpath(DA, "metric.tsv")

const J2 = joinpath(DA, "set.json")

const D1, D2, D3, D4, D5 =
    (mkpath(joinpath(TE, ba)) for ba in ("da", "us", "me", "us2", "me2"))

const T4 = joinpath(D3, "metric.tsv")

function rea(di)

    Nucleus.Table.rea(joinpath(di, "result.tsv"))

end

# ---- #
# TODO

GSEA.data_rank(D1, T2, J2; minimum = 15, maximum = 500)

const A1 = rea(D1)

@test size(A1) === (UM, 10)

A1 = A1[findall(nu_ -> !any(isnan, nu_), eachrow(Matrix(A1[!, 2:end]))), :]

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

GSEA.user_rank(D2, T3, J2)

const A2 = rea(D2)

@test size(A2) === (UM, 5)

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

    @test is_egal(collect(A2[id, 2:end]), r2)

end

# ---- #
# TODO

GSEA.metric_rank(D3, T1, T2, J2)

const A3 = Nucleus.Table.rea(T4)

@test size(A3) === (1000, 2)

@test is_egal(names(A3), ["Feature", "signal-to-noise-ratio"])

@test is_egal(sort(A3, 2)[[1, end], 2], [-1.8372355409610066, 1.7411005104346835])

const A4 = rea(D3)

@test size(A4) === (UM, 5)

# ---- #

GSEA.user_rank(D4, T4, J2)

const A5 = rea(D4)

@test is_egal(A5[!, 1:2], A4[!, 1:2])

# ---- #

GSEA.metric_rank(D5, T1, T2, J2; permutation = "set")

@test is_egal(rea(D5), A5)
