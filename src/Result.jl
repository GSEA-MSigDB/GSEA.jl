module Result

using Nucleus

using ..GSEA

function writ(di, al, s1_, nu_, s2_, st__, en_, R, um, s3_, t1, t2)

    i1_ = findall(!isnan, en_)

    s2_ = s2_[i1_]

    st__ = st__[i1_]

    en_ = en_[i1_]

    R = R[i1_, :]

    E = Matrix{Float64}(undef, lastindex(i1_), 4)

    E[:, 1] = en_

    GSEA.Normalization.make!(al, en_, R)

    E[:, 2] = en_

    i2_, pv_, qv_ = Nucleus.Significance.make(en_, R)

    E[i2_, 3] = pv_

    E[i2_, 4] = qv_

    Nucleus.Table.writ(
        joinpath(di, "result.tsv"),
        Nucleus.Table.make(
            "Set",
            s2_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
            E,
        ),
    )

    for id in unique!(
        vcat(Nucleus.Extreme.index(E[:, 1], um), filter!(!isnothing, indexin(s3_, s2_))),
    )

        st = s2_[id]

        GSEA.Plot.writ(
            joinpath(di, "$(Nucleus.Numbe.text_2(E[id, 1])).$st.html"),
            al,
            s1_,
            nu_,
            st__[id],
            Dict("title" => Dict("text" => st));
            t1,
            t2,
        )

    end

end

end
