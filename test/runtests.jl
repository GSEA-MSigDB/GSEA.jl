using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

for na in (
    "Algorithm",
    "Enrichment",
    "File",
    "Interface",
    "Normalization",
    "Plot",
    "Rando",
    "Result",
    "Sort",
)

    @info "🎬 Testing $na"

    run(`julia --project $na.jl`)

end

@info "🎬 Testing GSEA"

# ---- #

using Nucleus

include("_.jl")

# ---- #

const T1 = joinpath(TE, "_.tsv")

# ---- #

for (ba, re) in (
    ("1.cls", (1, 7)),
    ("GSE76137.cls", (1, 7)),
    ("CCLE_mRNA_20Q2_no_haem_phen.cls", (1, 900)),
)

    GSEA.cls(T1, joinpath(DA, ba))

    @test size(Nucleus.Table.rea(T1)) === re

end

# ---- #

for (ba, re) in (("1.gct", (13321, 190)),)

    GSEA.gct(T1, joinpath(DA, ba))

    @test size(Nucleus.Table.rea(T1)) === re

end

# ---- #

const J1 = joinpath(TE, "_.json")

for (ba, re) in (("1.gmt", 50), ("2.gmt", 5529))

    GSEA.gmt(J1, joinpath(DA, ba))

    @test length(Nucleus.Dictionary.rea(J1)) === re

end

# ---- #

const I1 = mkpath(joinpath(TE, "data_rank"))

const I2 = mkpath(joinpath(TE, "user_rank"))

const I3 = mkpath(joinpath(TE, "metric_rank"))

const J2 = joinpath(DA, "set.json")

const T2 = joinpath(DA, "data.tsv")

# ---- #

GSEA.data_rank(I1, T2, J2; minimum = 15, maximum = 500)

const A1 = Nucleus.Table.rea(joinpath(I1, "result.tsv"))

@test A1[!, 1] == [
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
]

@test findmax(Matrix(A1[!, 2:end])) === (0.756249206577638, CartesianIndex(2, 2))

# ---- #

function test(an, um)

    @test size(an, 1) === um

end

# ---- #

GSEA.user_rank(
    I2,
    joinpath(DA, "metric.tsv"),
    J2;
    more_plots = "HALLMARK_MYC_TARGETS_V1;HALLMARK_UV_RESPONSE_DN;HALLMARK_UV_RESPONSE_UP;ALIEN",
)

const A2 = Nucleus.Table.rea(joinpath(I2, "result.tsv"))

test(A2, 50)

for (id, r1, r2, r3, r4) in (
    (45, "HALLMARK_PANCREAS_BETA_CELLS", -0.35266, -1.36616, 0.0200837),
    (33, "HALLMARK_PROTEIN_SECRETION", -0.272096, -1.25207, 0.0686192),
    (36, "HALLMARK_MYC_TARGETS_V1", 0.603356, 2.73998, 0.000262812),
    (10, "HALLMARK_MYC_TARGETS_V2", 0.866579, 3.36557, 0.000262812),
)

    @test A2[id, 1] === r1

    @test isapprox(A2[id, 2], r2; atol = 1e-6)

    # TODO: Investigate.

    #@test isapprox(A2[id, 3], r3; atol = 1e-5)

    #@test isapprox(A2[id, 4], r4; atol = 1e-7)

end

# ---- #

const KE_ = (minimum = 15, maximum = 500, more_plots = "HALLMARK_UV_RESPONSE_DN;")

# ---- #

GSEA.metric_rank(I3, joinpath(DA, "target.tsv"), T2, J2; KE_...)

const A3 = Nucleus.Table.rea(joinpath(I3, "metric.tsv"))

@test size(A3) === (1000, 2)

@test names(A3) == ["Feature", "signal-to-noise-ratio"]

@test isapprox(sort(A3, 2)[[1, end], 2], [-1.8372355409610066, 1.7411005104346835])

const A4 = Nucleus.Table.rea(joinpath(I3, "result.tsv"))

test(A4, 8)

# ---- #

GSEA.user_rank(I2, joinpath(I3, "metric.tsv"), J2; KE_...)

const A5 = Nucleus.Table.rea(joinpath(I2, "result.tsv"))

@test A5[!, 1:2] == A4[!, 1:2]

# ---- #

GSEA.metric_rank(I3, joinpath(DA, "target.tsv"), T2, J2; permutation = "set", KE_...)

const A6 = Nucleus.Table.rea(joinpath(I3, "result.tsv"))

@test A6 == A5
