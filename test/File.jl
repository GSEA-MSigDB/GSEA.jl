using Test: @test

using GSEA

include("_.jl")

# ---- #

# 10.500 μs (114 allocations: 6.84 KiB)
# 10.708 μs (116 allocations: 7.17 KiB)
# 271.958 μs (7160 allocations: 497.98 KiB)

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

    A = GSEA.File.read_cls(cl)

    @btime GSEA.File.read_cls($cl)

    @test is_egal(names(A), vcat("Phenotype", map(id -> "Sample $id", 1:(size(A, 2) - 1))))

    @test A[!, 1][] === r1

    @test is_egal(Vector(A[1, 2:(1 + lastindex(r2))]), r2)

end

# ---- #

for (ba, re) in (("1.gct", (13321, 190)),)

    @test size(GSEA.File.read_gct(joinpath(DA, ba))) === re

end

# ---- #

# 181.959 μs (7984 allocations: 1.16 MiB)
# 15.765 ms (537839 allocations: 66.19 MiB)

for (ba, re) in (("1.gmt", 50), ("2.gmt", 5529))

    gm = joinpath(DA, ba)

    @test length(GSEA.File.read_gmt(gm)) === re

    @btime GSEA.File.read_gmt($gm)

end
