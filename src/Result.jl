module Result

using Public

using ..GSEA

function writ(di, al, s1_, nu_, s2_, st__, en_, R, um, s3_, t1, t2, t3)

    N = Matrix{Float64}(undef, lastindex(s2_), 4)

    N[:, 1] = en_

    GSEA.Normalization.make!(en_, R)

    N[:, 2] = en_

    in_, pv_, qv_ = Public.number_significance(en_, R)

    N[in_, 3] = pv_

    N[in_, 4] = qv_

    Public.write_table(
        joinpath(di, "result.tsv"),
        Public.make_table(
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

    for nd in unique!(
        vcat(
            Public.index_extreme(N[:, 1], um),
            filter!(!isnothing, indexin(s3_, s2_)),
        ),
    )

        st = s2_[nd]

        GSEA.Plot.writ(
            joinpath(di, "$(Public.text_2(N[nd, 1])).$st.html"),
            al,
            s1_,
            nu_,
            st__[nd],
            Dict("title" => Dict("text" => st));
            t1,
            t2,
            t3,
        )

    end

end

end
