using Test: @test

include("_.jl")

# ---- #

function test(fu, ba, re)

    ts = joinpath(TE, "_.tsv")

    fu(ts, joinpath(DA, ba))

    @test size(Nucleus.Table.rea(ts)) === re

end

# ---- #

for (ba, re) in (
    ("1.cls", (1, 7)),
    ("GSE76137.cls", (1, 7)),
    ("CCLE_mRNA_20Q2_no_haem_phen.cls", (1, 900)),
)

    test(GSEA.cls, ba, re)

end

# ---- #

for (ba, re) in (("1.gct", (13321, 190)),)

    test(GSEA.gct, ba, re)

end

# ---- #

const JS = joinpath(TE, "_.json")

for (ba, re) in (("h.all.v7.1.symbols.gmt", 50), ("c2.all.v7.1.symbols.gmt", 5529))

    GSEA.gmt(JS, joinpath(DA, ba))

    @test length(Nucleus.Dictionary.rea(JS)) === re

end
