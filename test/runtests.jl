using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Nucleus

# ---- #

const DI = pkgdir(GSEA, "data")

# ---- #

const TE = joinpath(tempdir(), "GSEA")

rm(TE; recursive = true, force = true)

mkdir(TE)

# ---- #

const TS = joinpath(TE, "_.tsv")

# ---- #

for (cl, re) in (
    ("1.cls", (1, 7)),
    ("GSE76137.cls", (1, 7)),
    ("CCLE_mRNA_20Q2_no_haem_phen.cls", (1, 900)),
)

    GSEA.cls(TS, joinpath(DI, cl))

    @test size(Nucleus.Table.rea(TS)) === re

end

# ---- #

for (gc, re) in (("1.gct", (13321, 190)),)

    GSEA.gct(TS, joinpath(DI, gc))

    @test size(Nucleus.Table.rea(TS)) === re

end

# ---- #

const JS = joinpath(TE, "_.json")

for (gm, re) in (("1.gmt", 50), ("2.gmt", 5529))

    GSEA.gmt(JS, joinpath(DI, gm))

    @test length(Nucleus.Dictionary.rea(JS)) === re

end

# ---- #

@test GSEA._select_sort('a':'f', [1, NaN, 3, NaN, 5]) == (['e', 'c', 'a'], [5.0, 3.0, 1.0])

# ---- #

const SE_, ME___ = GSEA._separat(GSEA.read_gmt(joinpath(DI, "h.all.v7.1.symbols.gmt")))

# ---- #

# 3.022 ms (108 allocations: 934.22 KiB)
# 2.649 ms (108 allocations: 934.22 KiB)
# 17.202 ms (108 allocations: 934.22 KiB)
# 17.186 ms (108 allocations: 934.22 KiB)
# 10.127 ms (108 allocations: 934.22 KiB)
# 9.001 ms (108 allocations: 934.22 KiB)
#
# 3.051 ms (32 allocations: 1.86 MiB)
# 2.790 ms (32 allocations: 1.86 MiB)
# 21.514 ms (32 allocations: 1.86 MiB)
# 21.662 ms (32 allocations: 1.86 MiB)
# 13.287 ms (32 allocations: 1.86 MiB)
# 11.667 ms (32 allocations: 1.86 MiB)

for al in AL_

    en_ = GSEA.enrich(al, A1_, N3_, ME___)

    @test !issorted(en_)

    se_, en_ = GSEA._select_sort(SE_, en_)

    @test se_[1:2] == ["HALLMARK_MYC_TARGETS_V2", "HALLMARK_MYC_TARGETS_V1"]

    #@btime GSEA.enrich($al, A1_, N3_, ME___)

end

# ---- #

const DG = joinpath(tempdir(), "GSEA")

if isdir(DG)

    rm(DG; recursive = true)

end

mkdir(DG)

const FS = joinpath(DI, "set.json")

const FD = joinpath(DI, "data.tsv")

# ---- #

for (al, re) in zip(AL_, ("KS", "KSa", "KLioM", "KLioP", "KLi", "KLi1"))

    @test GSEA._set_algorithm(lowercase(re)) === al

end

# ---- #

const OD = mkpath(joinpath(DG, "data_rank"))

const SD_, ED = GSEA.data_rank(OD, FD, FS; minimum_set_size = 15, maximum_set_size = 500)

@test isfile(joinpath(OD, "enrichment.tsv"))

@test SD_ == [
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_KRAS_SIGNALING_DN",
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
]

@test findmax(ED) === (0.756249206577638, CartesianIndex(2, 2))

# ---- #

# 435.333 μs (1802 allocations: 3.03 MiB)
# 463.708 μs (1502 allocations: 3.69 MiB)

for al in (AL_[1], AL_[end])

    seed!(20231103)

    en_ = randn(100)

    ra = randn(100, 1000)

    GSEA._normalize_enrichment!(al, en_, ra)

    #@btime GSEA._normalize_enrichment!($al, $en_, $ra)

end

# ---- #

function test_result(an, US)

    @test size(an, 1) === US

end

# ---- #

const OU = mkpath(joinpath(DG, "user_rank"))

GSEA.user_rank(
    OU,
    joinpath(DI, "metric.tsv"),
    FS;
    more_sets_to_plot = "HALLMARK_MYC_TARGETS_V1;HALLMARK_UV_RESPONSE_DN;HALLMARK_UV_RESPONSE_UP;ALIEN",
)

const RU = Nucleus.Table.rea(joinpath(OU, "result.tsv"))

test_result(RU, 50)

for (id, r1, r2, r3, r4) in (
    (45, "HALLMARK_PANCREAS_BETA_CELLS", -0.35266, -1.36616, 0.0200837),
    (33, "HALLMARK_PROTEIN_SECRETION", -0.272096, -1.25207, 0.0686192),
    (36, "HALLMARK_MYC_TARGETS_V1", 0.603356, 2.73998, 0.000262812),
    (10, "HALLMARK_MYC_TARGETS_V2", 0.866579, 3.36557, 0.000262812),
)

    @test RU[id, 1] === r1

    @test isapprox(RU[id, 2], r2; atol = 1e-6)

    # TODO: Investigate.

    #@test isapprox(RU[id, 3], r3; atol = 1e-5)

    #@test isapprox(RU[id, 4], r4; atol = 1e-7)

end

# ---- #

const OM = mkpath(joinpath(DG, "metric_rank"))

GSEA.metric_rank(
    OM,
    joinpath(DI, "target.tsv"),
    FD,
    FS;
    minimum_set_size = 15,
    maximum_set_size = 500,
)

const ME = Nucleus.Table.rea(joinpath(OM, "metric.tsv"))

@test size(ME) === (1000, 2)

@test names(ME) == ["Feature", "signal-to-noise-ratio"]

@test isapprox(sort(ME, 2)[[1, end], 2], [-1.8372355409610066, 1.7411005104346835])

const RU = Nucleus.Table.rea(joinpath(OM, "result.tsv"))

test_result(RU, 8)
