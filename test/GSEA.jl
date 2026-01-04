using Test: @test

include("_.jl")

# ---- #

function test(fu, ba, re)

    ts = joinpath(TE, "_.tsv")

    fu(ts, joinpath(DA, ba))

    @test size(Public.make_part(Public.read_table(ts))[4]) === re

end

# ---- #

for (ba, re) in (
    ("1.cls", (1, 6)),
    ("GSE76137.cls", (1, 6)),
    ("CCLE_mRNA_20Q2_no_haem_phen.cls", (1, 899)),
)

    test(GSEA.cls, ba, re)

end

# ---- #

for (ba, re) in (("1.gct", (13321, 189)),)

    test(GSEA.gct, ba, re)

end

# ---- #

const UM = 50

# ---- #

const J1 = joinpath(TE, "_.json")

for (ba, re) in
    (("h.all.v7.1.symbols.gmt", UM), ("c2.all.v7.1.symbols.gmt", 5529))

    GSEA.gmt(J1, joinpath(DA, ba))

    @test length(Public.read_pair(J1)) === re

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

    Public.make_part(Public.read_table(joinpath(di, "result.tsv")))

end

# ---- #

GSEA.data_rank(D1, T2, J2; minimum = 15, maximum = 500)

const _, S6_, _, N1 = rea(D1)

@test size(N1) === (UM, 9)

const IN_ = findall(nu_ -> !any(isnan, nu_), eachrow(N1))

@test is_egal(
    S6_[IN_],
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

@test findmax(N1[IN_, :]) === (0.756249206577638, CartesianIndex(2, 2))

# ---- #

GSEA.user_rank(D2, T3, J2)

const _, S7_, _, N2 = rea(D2)

@test size(N2) === (UM, 4)

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

    @test S7_[nd] === r1

    @test N2[nd, :] == r2

end

# ---- #

GSEA.metric_rank(D3, T1, T2, J2)

const ST, _, S8_, N3 = Public.make_part(Public.read_table(T4))

@test ST === "Gene"

@test S8_[] === "signal-to-noise-ratio"

const N3_ = sort!(vec(N3))

@test lastindex(N3_) === 1000

@test N3_[1] === -1.8372355409610066

@test N3_[end] === 1.7411005104346835

const _, S9_, _, N4 = rea(D3)

@test size(N4) === (UM, 4)

# ---- #

GSEA.user_rank(D4, T4, J2)

const _, S0_, _, N5 = rea(D4)

@test is_egal(S9_, S0_)

@test is_egal(N4[:, 1], N5[:, 1])

# ---- #

GSEA.metric_rank(D5, T1, T2, J2; permutation = "set")

const _, S9_, _, N4 = rea(D4)

@test is_egal(S9_, S0_)

@test is_egal(N4, N5)
