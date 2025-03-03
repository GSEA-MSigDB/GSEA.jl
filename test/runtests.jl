using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Random: seed!

using Nucleus

# ---- #

const DI = pkgdir(GSEA, "data")

const TE = joinpath(tempdir(), "GSEA")

rm(TE; recursive = true, force = true)

mkdir(TE)

# ---- #

# 10.000 μs (114 allocations: 6.84 KiB)
# 10.291 μs (116 allocations: 7.17 KiB)
# 403.875 μs (7160 allocations: 473.16 KiB)

for (cl, na, re) in (
    ("1.cls", "CNTRL_LPS", [1, 1, 1, 2, 2, 2]),
    ("GSE76137.cls", "Proliferating_Arrested", [1, 2, 1, 2, 1, 2]),
    (
        "CCLE_mRNA_20Q2_no_haem_phen.cls",
        "HER2",
        [1.087973, -1.358492, -1.178614, -0.77898, 0.157222, 1.168224, -0.360195, 0.608629],
    ),
)

    cl = joinpath(DI, cl)

    an = GSEA.read_cls(cl)

    @btime GSEA.read_cls($cl)

    na_ = names(an)

    N = Matrix(an[:, 2:end])

    @test na_[1] === "Phenotype"

    @test an[:, 1][] === na

    @test all(startswith("Sample "), na_[2:end])

    @test eltype(N) === eltype(re)

    @test N[1, eachindex(re)] == re

end

# ---- #

# 100.965 ms (71705 allocations: 23.67 MiB)

for (gc, re) in (("1.gct", (13321, 190)),)

    gc = joinpath(DI, gc)

    @test size(GSEA.read_gct(gc)) === re

    @btime GSEA.read_gct($gc)

end

# ---- #

# 287.417 μs (7984 allocations: 1.12 MiB)
# 21.892 ms (537839 allocations: 62.61 MiB)

for (gm, re) in (("1.gmt", 50), ("2.gmt", 5529))

    gm = joinpath(DI, gm)

    di = GSEA.read_gmt(gm)

    @btime GSEA.read_gmt($gm)

    @test typeof(di) === Dict{String, Vector{String}}

    @test length(di) === re

end

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

const AL_ = GSEA.KS(), GSEA.KSa(), GSEA.KLioM(), GSEA.KLioP(), GSEA.KLi(), GSEA.KLi1()

# ---- #

for (al, re) in zip(AL_, ("KS", "KSa", "KLioM", "KLioP", "KLi", "KLi1"))

    @test GSEA.text(al) === re

end

# ---- #

const N1_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const B1_ = [true, false, true, false, true, true, false, false, true]

const N2_ = randn(100000)

const B2_ = rand(Bool, lastindex(N2_))

# ---- #

# 70.385 ns (0 allocations: 0 bytes)
# 70.385 ns (0 allocations: 0 bytes)
# 8.291 ns (0 allocations: 0 bytes)
# 19.767 ns (0 allocations: 0 bytes)
# 309.958 μs (0 allocations: 0 bytes)

const N0 = -0.25

for (nu_, ex, bo_, re) in (
    (N1_, 0.1, B1_, (N0, 0.24581982412836917)),
    (N1_, 0.5, B1_, (N0, 0.21402570288861142)),
    (N1_, 1, B1_, (N0, 0.15625)),
    (N1_, 2, B1_, (N0, 0.06226650062266501)),
    (N2_, 1, B2_, nothing),
)

    al = GSEA.KS()

    @test isnothing(re) || GSEA.make_normalizer(al, nu_, ex, bo_) === re

    @btime GSEA.make_normalizer($al, $nu_, $ex, $bo_)

end

# ---- #

# 117.039 ns (0 allocations: 0 bytes)
# 117.011 ns (0 allocations: 0 bytes)
# 7.375 ns (0 allocations: 0 bytes)
# 25.016 ns (0 allocations: 0 bytes)
# 92.875 μs (0 allocations: 0 bytes)

for (nu_, ex, bo_, re) in (
    (N1_, 0.1, B1_, (0.24581982412836917, 0.14006007078470165)),
    (N1_, 0.5, B1_, (0.21402570288861142, 0.12366213677204271)),
    (N1_, 1, B1_, (0.15625, 0.09615384615384615)),
    (N1_, 2, B1_, (0.06226650062266501, 0.04533091568449683)),
    (N2_, 1, B2_, nothing),
)

    al = GSEA.KLioM()

    @test isnothing(re) || GSEA.make_normalizer(al, nu_, ex, bo_) === re

    @btime GSEA.make_normalizer($al, $nu_, $ex, $bo_)

end

# ---- #

for (n1, n2, re) in ((1 / 3, 0.5, -1.0),)

    @test GSEA.make_normalizer(n1, n2) === re

end

# ---- #

const A1_, N3_ =
    eachcol(reverse!(Nucleus.Table.rea(joinpath(DI, "myc.tsv"); select = [1, 2])))

const M1_ = GSEA.read_gmt(joinpath(DI, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]

# ---- #

# 19.433 ns (0 allocations: 0 bytes)
# 17.869 ns (0 allocations: 0 bytes)
# 225.103 ns (0 allocations: 0 bytes)
# 225.069 ns (0 allocations: 0 bytes)
# 126.633 ns (0 allocations: 0 bytes)
# 118.359 ns (0 allocations: 0 bytes)
# 45.208 μs (0 allocations: 0 bytes)
# 37.500 μs (0 allocations: 0 bytes)
# 325.833 μs (0 allocations: 0 bytes)
# 326.042 μs (0 allocations: 0 bytes)
# 186.208 μs (0 allocations: 0 bytes)
# 164.833 μs (0 allocations: 0 bytes)
#
# 16.825 ns (0 allocations: 0 bytes)
# 16.325 ns (0 allocations: 0 bytes)
# 281.034 ns (0 allocations: 0 bytes)
# 281.034 ns (0 allocations: 0 bytes)
# 155.565 ns (0 allocations: 0 bytes)
# 155.050 ns (0 allocations: 0 bytes)
# 45.833 μs (0 allocations: 0 bytes)
# 37.541 μs (0 allocations: 0 bytes)
# 410.750 μs (0 allocations: 0 bytes)
# 410.833 μs (0 allocations: 0 bytes)
# 243.375 μs (0 allocations: 0 bytes)
# 222.541 μs (0 allocations: 0 bytes)

for (na_, nu_, me_, re_) in (
    (
        ['K', 'Q', 'J', 'X', '9', '8', '7', '6', '5', '4', '3', '2', 'A'],
        [6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6.0],
        ['K', 'A'],
        (-0.5, 0.0, 0.0, 0.0, 0.0, 0.0),
    ),
    (
        A1_,
        N3_,
        M1_,
        (
            0.7651927829281453,
            0.41482514169516305,
            1.1181841586127337,
            1.1140922794954267,
            1.1161382190540838,
            1.2297916337424049,
        ),
    ),
)

    bo_ = Nucleus.Collection.is_in(na_, me_)

    for (al, re) in zip(AL_, re_)

        ex = 1

        @test isapprox(GSEA._enrich!(al, nu_, ex, bo_, nothing), re; atol = 1e-11)

        @btime GSEA._enrich!($al, $nu_, $ex, $bo_, nothing)

        GSEA.plot("", al, na_, nu_, me_, Dict("title" => Dict("text" => GSEA.text(al)));)

    end

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
