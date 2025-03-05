using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Nucleus

include("_.jl")

# ---- #

# 10.542 μs (114 allocations: 6.84 KiB)
# 10.875 μs (116 allocations: 7.17 KiB)
# 412.208 μs (7160 allocations: 473.16 KiB)

for (cl, r1, r2) in (
    ("1.cls", "CNTRL_LPS", [1, 1, 1, 2, 2, 2]),
    ("GSE76137.cls", "Proliferating_Arrested", [1, 2, 1, 2, 1, 2]),
    (
        "CCLE_mRNA_20Q2_no_haem_phen.cls",
        "HER2",
        [1.087973, -1.358492, -1.178614, -0.77898, 0.157222, 1.168224, -0.360195, 0.608629],
    ),
)

    cl = joinpath(DI, cl)

    an = GSEA.File.read_cls(cl)

    #@btime GSEA.File.read_cls($cl)

    @test names(an) == vcat("Phenotype", map(id -> "Sample $id", 1:(size(an, 2) - 1)))

    @test an[!, 1][] === r1

    N = Matrix(an[!, 2:end])

    @test eltype(N) === eltype(r2)

    @test N[1, eachindex(r2)] == r2

end

# ---- #

# 103.396 ms (71705 allocations: 23.67 MiB)

for (gc, re) in (("1.gct", (13321, 190)),)

    gc = joinpath(DI, gc)

    @test size(GSEA.File.read_gct(gc)) === re

    #@btime GSEA.File.read_gct($gc)

end

# ---- #

# 285.959 μs (7984 allocations: 1.12 MiB)
# 22.009 ms (537839 allocations: 62.61 MiB)

for (gm, re) in (("1.gmt", 50), ("2.gmt", 5529))

    gm = joinpath(DI, gm)

    di = GSEA.File.read_gmt(gm)

    #@btime GSEA.File.read_gmt($gm)

    @test length(di) === re

end
