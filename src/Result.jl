module Result

using Nucleus

using ..GSEA

function writ(di, al, n1_, nu_, n3_, n2__, en_, R, um, n4_, a1, a2)

    in_ = findall(!isnan, en_)

    n3_ = n3_[in_]

    n2__ = n2__[in_]

    en_ = en_[in_]

    R = R[in_, :]

    E = Matrix{Float64}(undef, lastindex(in_), 4)

    E[:, 1] = en_

    GSEA.Normalization.make!(al, en_, R)

    E[:, 2] = en_

    in_, pv_, qv_ = Nucleus.Significance.make(en_, R)

    E[in_, 3] = pv_

    E[in_, 4] = qv_

    Nucleus.Table.writ(
        joinpath(di, "result.tsv"),
        Nucleus.Table.make(
            "Set",
            n3_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
            E,
        ),
    )

    for id in unique!(
        vcat(Nucleus.Extreme.index(E[:, 1], um), filter!(!isnothing, indexin(n4_, n3_))),
    )

        n3 = n3_[id]

        GSEA.Plot.writ(
            joinpath(di, "$(Nucleus.Numbe.text_2(E[id, 1])).$n3.html"),
            al,
            n1_,
            nu_,
            n2__[id],
            Dict("title" => Dict("text" => n3));
            a1,
            a2,
        )

    end

end

end
