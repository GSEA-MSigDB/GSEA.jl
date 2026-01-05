using Public

using GSEA

########################################

const AL_ = GSEA.KS0(), GSEA.A0(), GSEA.DA2(), GSEA.DA2W(), GSEA.DA2W0W()

########################################

const S1_, N1_ = GSEA.make_sort(
    [
        "Aa",
        "22",
        "33",
        "44",
        "55",
        "66",
        "77",
        "88",
        "99",
        "Xx",
        "Jj",
        "Qq",
        "Kk",
    ],
    [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6.0],
)

for st_ in (["Aa", "Kk", "__"], ["66", "77", "88", "__"]), al in AL_

    GSEA.write_enrichment(
        "",
        al,
        S1_,
        N1_,
        st_,
        Dict("title" => Dict("text" => "$al")),
    )

end

########################################

const PA = cp(pkgdir(GSEA, "ou"), joinpath(tempdir(), "GSEA"); force = true)

const S2_, N2_ = GSEA.make_sort(
    eachcol(Public.read_table(joinpath(GSEA.P1, "myc.tsv"))[!, 1:2])...,
)

const N = hcat(N2_, N2_ * 2)

const DI = GSEA.read_gmt(joinpath(GSEA.P1, "h.all.v7.1.symbols.gmt"))

const S3_ = collect(keys(DI))

const ST__ = collect(values(DI))

for al in AL_

    st = "$al"

    GSEA.write_enrichment(
        joinpath(PA, st),
        al,
        S2_,
        ("Sa1", "Sa2"),
        N,
        S3_,
        ST__,
        reduce(
            hcat,
            GSEA.number_enrichment(al, S2_, nu_, ST__) for nu_ in eachcol(N)
        ),
        Dict("title" => Dict("text" => st));
        um = 1,
        t3 = "High expression",
    )

end
