using Test: @test

include("_.jl")

# ---- #

# 10.959 μs (114 allocations: 6.84 KiB)
# 10.875 μs (116 allocations: 7.17 KiB)
# 279.417 μs (7160 allocations: 497.98 KiB)

for (ba, r1, r2) in (
    ("1.cls", "CNTRL_LPS", [1, 1, 1, 2, 2, 2]),
    ("GSE76137.cls", "Proliferating_Arrested", [1, 2, 1, 2, 1, 2]),
    (
        "CCLE_mRNA_20Q2_no_haem_phen.cls",
        "HER2",
        [1.087973, -1.358492, -1.178614, -0.77898, 0.157222, 1.168224, -0.360195, 0.608629],
    ),
)

    cl = joinpath(DA, ba)

    st, s1_, s2_, N = Nucleus.Table.ge(GSEA.File.read_cls(cl))

    #@btime GSEA.File.read_cls($cl)

    @test st === "Phenotype"

    @test s1_[] === r1

    @test is_egal(s2_, map(id -> "Sample $id", eachindex(s2_)))

    @test is_egal(N[eachindex(r2)], r2)

end

# ---- #

for (ba, re) in (("1.gct", (13321, 190)),)

    @test size(GSEA.File.read_gct(joinpath(DA, ba))) === re

end

# ---- #

# 182.166 μs (7986 allocations: 1.16 MiB)
# 14.968 ms (537841 allocations: 66.19 MiB)

for (ba, re) in (("h.all.v7.1.symbols.gmt", 50), ("c2.all.v7.1.symbols.gmt", 5529))

    gm = joinpath(DA, ba)

    @test length(GSEA.File.read_gmt(gm)) === re

    #@btime GSEA.File.read_gmt($gm)

end
