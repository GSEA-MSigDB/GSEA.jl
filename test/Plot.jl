using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Nucleus

include("_.jl")

# ---- #

for (n1_, nu_, n2_) in (
        (
            ['K', 'Q', 'J', 'X', '9', '8', '7', '6', '5', '4', '3', '2', 'A'],
            [6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6.0],
            ['K', 'A'],
        ),
        (
            eachcol(
                reverse!(Nucleus.Table.rea(joinpath(DI, "myc.tsv"); select = [1, 2])),
            )...,
            GSEA.File.read_gmt(joinpath(DI, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"],
        ),
    ),
    al in AL_

    GSEA.Plot.writ("", al, n1_, nu_, n2_, Dict("title" => Dict("text" => string(al))))

end
