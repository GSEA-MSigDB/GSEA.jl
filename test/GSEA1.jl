using Test: @test

using Nucleus

using GSEA

include("_.jl")

# ---- #

const TS = joinpath(TE, "_.tsv")

# ---- #

for (ba, re) in (
    ("1.cls", (1, 7)),
    ("GSE76137.cls", (1, 7)),
    ("CCLE_mRNA_20Q2_no_haem_phen.cls", (1, 900)),
)

    GSEA.cls(TS, joinpath(DA, ba))

    @test size(Nucleus.Table.rea(TS)) === re

end

# ---- #

for (ba, re) in (("1.gct", (13321, 190)),)

    GSEA.gct(TS, joinpath(DA, ba))

    @test size(Nucleus.Table.rea(TS)) === re

end

# ---- #

const JS = joinpath(TE, "_.json")

for (ba, re) in (("1.gmt", 50), ("2.gmt", 5529))

    GSEA.gmt(JS, joinpath(DA, ba))

    @test length(Nucleus.Dictionary.rea(JS)) === re

end
