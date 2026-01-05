using Test: @test

using GSEA

########################################

for (s1, r1, r2) in (
    ("1.cls", "CNTRL vs LPS", [1, 1, 1, 2, 2, 2]),
    ("GSE76137.cls", "Proliferating vs Arrested", [1, 2, 1, 2, 1, 2]),
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

    _, s2, _, N = GSEA.read_cls(joinpath(GSEA.P1, s1))

    @test s2 === r1

    @test N[eachindex(r2)] == r2

end

########################################

@test size(GSEA.read_gct(joinpath(GSEA.P1, "1.gct"))[4]) === (13321, 189)

########################################

for (st, re) in
    (("h.all.v7.1.symbols.gmt", 50), ("c2.all.v7.1.symbols.gmt", 5529))

    @test length(GSEA.read_gmt(joinpath(GSEA.P1, st))) === re

end
