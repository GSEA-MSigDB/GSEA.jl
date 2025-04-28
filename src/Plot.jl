module Plot

using Nucleus

using ..GSEA

function writ(ht, al, s1_, nu_, s2_, la = Dict{String, Any}(); t1 = "Low", t2 = "High")

    s1_, nu_ = GSEA.Sort.make(s1_, nu_)

    um = lastindex(s1_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    xc_ = 1:um

    i1_ = findall(<(0), nu_)

    i2_ = findall(>=(0), nu_)

    bo_ = Nucleus.Collection.is_in(s1_, s2_)

    cu_ = Vector{Float64}(undef, um)

    en = GSEA.Enrichment.make!(al, nu_, bo_, cu_)

    an = Dict(
        "y" => 0,
        "font" => Dict("size" => 16),
        "borderpad" => 4.8,
        "borderwidth" => 2.64,
        "bordercolor" => Nucleus.Color.LI,
        "showarrow" => false,
    )

    po = um * 0.008

    Nucleus.Plotly.writ(
        ht,
        (
            merge(
                tr,
                Dict(
                    "y" => nu_[i1_],
                    "x" => xc_[i1_],
                    "text" => s1_[i1_],
                    "fillcolor" => Nucleus.Color.BL,
                ),
            ),
            merge(
                tr,
                Dict(
                    "y" => nu_[i2_],
                    "x" => xc_[i2_],
                    "text" => s1_[i2_],
                    "fillcolor" => Nucleus.Color.RE,
                ),
            ),
            Dict(
                "yaxis" => "y2",
                "y" => zeros(sum(bo_)),
                "x" => xc_[bo_],
                "mode" => "markers",
                "marker" => Dict(
                    "symbol" => "line-ns",
                    "size" => 24,
                    "line" => Dict("width" => 2, "color" => "#000000cc"),
                ),
                "text" => s1_[bo_],
                "hoverinfo" => "x+text",
            ),
            merge(
                tr,
                Dict(
                    "yaxis" => "y3",
                    "y" => cu_,
                    "x" => xc_,
                    "text" => s1_,
                    "fillcolor" => "#07fa07",
                ),
            ),
        ),
        Nucleus.Dictionary.make(
            Dict(
                "showlegend" => false,
                "yaxis" => Dict("domain" => (0, 0.24), "title" => Dict("text" => "Score")),
                "yaxis2" => Dict(
                    "domain" => (0.248, 0.32),
                    "title" => Dict("text" => "Set"),
                    "tickvals" => (),
                ),
                "yaxis3" =>
                    Dict("domain" => (0.328, 1), "title" => Dict("text" => "Δ Enrichment")),
                "xaxis" => Dict(
                    "title" => Dict("text" => "Feature ($um)"),
                    "showspikes" => true,
                    "spikemode" => "across",
                    "spikedash" => "solid",
                    "spikethickness" => -1,
                    "spikecolor" => "#000000",
                ),
                "annotations" => (
                    Dict(
                        "yref" => "paper",
                        "xref" => "paper",
                        "y" => 1.056,
                        "text" => "Enrichment = <b>$(Nucleus.Numbe.text_4(en))</b>",
                        "font" => Dict("size" => 24, "color" => "#000000"),
                        "showarrow" => false,
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => 1 - po,
                            "xanchor" => "right",
                            "text" => t2,
                            "font" => Dict("color" => Nucleus.Color.RE),
                        ),
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => um + po,
                            "xanchor" => "left",
                            "text" => t1,
                            "font" => Dict("color" => Nucleus.Color.BL),
                        ),
                    ),
                ),
            ),
            la,
        ),
    )

end

function writ(fi, al, s2_, s1_, N, s3_, st__, E, um = 2, la = Dict{String, Any}(); ke_...)

    in_ = findall(en_ -> all(!isnan, en_), eachrow(E))

    s3_ = s3_[in_]

    st__ = st__[in_]

    E = E[in_, :]

    Nucleus.Table.writ("$fi.tsv", Nucleus.Table.make("Set", s3_, s2_, E))

    Nucleus.HeatPlot.writ(
        "$fi.html",
        s3_,
        s2_,
        E,
        Nucleus.Dictionary.make(
            Dict(
                "title" => Dict("text" => "Enrichment"),
                "yaxis" => Dict("title" => Dict("text" => "Set")),
                "xaxis" => Dict("title" => Dict("text" => "Sample")),
            ),
            la,
        ),
    )

    for in_ in CartesianIndices(E)[Nucleus.Extreme.index(vec(E), um)]

        i3, i2 = Tuple(in_)

        st = s3_[i3]

        writ(
            "$fi.$(Nucleus.Numbe.text_2(E[in_])).$(s2_[i2]).$st.html",
            al,
            s1_,
            N[:, i2],
            st__[i3],
            Dict("title" => Dict("text" => st));
            ke_...,
        )

    end

end

end
