using Test: @test

include("_.jl")

# ---- #

# 9.292 μs (78 allocations: 4.86 KiB)
# 9.666 μs (80 allocations: 5.19 KiB)
# 111.583 μs (4924 allocations: 283.03 KiB)

for (ba, r1, r2) in (
    ("1.cls", "CNTRL_LPS", [1, 1, 1, 2, 2, 2]),
    ("GSE76137.cls", "Proliferating_Arrested", [1, 2, 1, 2, 1, 2]),
    (
        "CCLE_mRNA_20Q2_no_haem_phen.cls",
        "HER2",
        [
            1.087973,
            -1.358492,
            -1.178614,
            -0.77898,
            0.157222,
            1.168224,
            -0.360195,
            0.608629,
        ],
    ),
)

    cl = joinpath(DA, ba)

    _, st, _, N = GSEA.File.read_cls(cl)

    #@btime GSEA.File.read_cls($cl)

    @test st === r1

    @test is_egal(N[eachindex(r2)], r2)

end

# ---- #

for (ba, re) in (("1.gct", (13321, 189)),)

    N = GSEA.File.read_gct(joinpath(DA, ba))[4]

    @test typeof(N) === Matrix{Float64}

    @test size(N) === re

end

# ---- #

# 184.125 μs (7986 allocations: 1.16 MiB)
# 15.507 ms (537841 allocations: 66.19 MiB)

for (ba, re) in
    (("h.all.v7.1.symbols.gmt", 50), ("c2.all.v7.1.symbols.gmt", 5529))

    gm = joinpath(DA, ba)

    @test length(GSEA.File.read_gmt(gm)) === re

    #@btime GSEA.File.read_gmt($gm)

end
