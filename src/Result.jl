module Result

using Nucleus

using ..GSEA

function writ(di, al, s1_, nu_, s2_, st__, en_, R, um, s3_, t1, t2, t3)

    N = Matrix{Float64}(undef, lastindex(s2_), 4)

    N[:, 1] = en_

    GSEA.Normalization.make!(en_, R)

    N[:, 2] = en_

    in_, pv_, qv_ = Nucleus.Significance.make(en_, R)

    N[in_, 3] = pv_

    N[in_, 4] = qv_

    Nucleus.Table.writ(
        joinpath(di, "result.tsv"),
        Nucleus.Table.make(
            "Set",
            s2_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
            N,
        ),
    )

    # TODO: Use normalized enrichment
    in_ = findall(!isnan, N[:, 1])

    s2_ = s2_[in_]

    st__ = st__[in_]

    N = N[in_, :]

    for id in unique!(
        vcat(Nucleus.Extreme.index(N[:, 1], um), filter!(!isnothing, indexin(s3_, s2_))),
    )

        st = s2_[id]

        GSEA.Plot.writ(
            joinpath(di, "$(Nucleus.Numbe.text_2(N[id, 1])).$st.html"),
            al,
            s1_,
            nu_,
            st__[id],
            Dict("title" => Dict("text" => st));
            t1,
            t2,
            t3,
        )

    end

end

end
