using Random: seed!

using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

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

    GSEA.CommandLineInterface.cls(T1, joinpath(DI, ba))

    @test size(Nucleus.Table.rea(T1)) === re

end

# ---- #

for (ba, re) in (("1.gct", (13321, 190)),)

    GSEA.CommandLineInterface.gct(T1, joinpath(DI, ba))

    @test size(Nucleus.Table.rea(T1)) === re

end

# ---- #

const J1 = joinpath(TE, "_.json")

for (ba, re) in (("1.gmt", 50), ("2.gmt", 5529))

    GSEA.CommandLineInterface.gmt(J1, joinpath(DI, ba))

    @test length(Nucleus.Dictionary.rea(J1)) === re

end

# ---- #

for (al, re) in zip(("ks", "ksa", "kliom", "kliop", "kli", "kli1"), AL_)

    @test GSEA.CommandLineInterface.make_algorithm(al) === re

end

# ---- #

const J2 = joinpath(DI, "set.json")

const T2 = joinpath(DI, "data.tsv")

# ---- #

const O1 = mkpath(joinpath(TE, "data_rank"))

GSEA.CommandLineInterface.data_rank(O1, T2, J2; minimum = 15, maximum = 500)

const A1 = Nucleus.Table.rea(joinpath(O1, "result.tsv"))

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

# 388.958 μs (1500 allocations: 2.26 MiB)
# 426.584 μs (1200 allocations: 2.91 MiB)

for al in (AL_[1], AL_[end])

    seed!(20231103)

    en_ = randn(100)

    ra = randn(100, 1000)

    GSEA.CommandLineInterface.make_normalized!(al, en_, ra)

    @btime GSEA.CommandLineInterface.make_normalized!($al, $en_, $ra)

end

# ---- #

function test_result(an, US)

    @test size(an, 1) === US

end

# ---- #

const O2 = mkpath(joinpath(TE, "user_rank"))

GSEA.CommandLineInterface.user_rank(
    O2,
    joinpath(DI, "metric.tsv"),
    J2;
    more_plots = "HALLMARK_MYC_TARGETS_V1;HALLMARK_UV_RESPONSE_DN;HALLMARK_UV_RESPONSE_UP;ALIEN",
)

const A2 = Nucleus.Table.rea(joinpath(O2, "result.tsv"))

test_result(A2, 50)

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

const O3 = mkpath(joinpath(TE, "metric_rank"))

GSEA.metric_rank(O3, joinpath(DI, "target.tsv"), T2, J2; minimum = 15, maximum = 500)

const A3 = Nucleus.Table.rea(joinpath(O3, "metric.tsv"))

@test size(A3) === (1000, 2)

@test names(A3) == ["Feature", "signal-to-noise-ratio"]

@test isapprox(sort(A3, 2)[[1, end], 2], [-1.8372355409610066, 1.7411005104346835])

const A4 = Nucleus.Table.rea(joinpath(O3, "result.tsv"))

test_result(A4, 8)
