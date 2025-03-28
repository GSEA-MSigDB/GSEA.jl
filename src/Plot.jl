module Plot

using Nucleus

using ..GSEA

function writ(ht, al, n1_, nu_, n2_, la = Dict{String, Any}(); a1 = "Low", a2 = "High")

    n1_, nu_ = GSEA.Sort.make(n1_, nu_)

    um = lastindex(n1_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    xc_ = 1:um

    i1_ = findall(<(0), nu_)

    i2_ = findall(>=(0), nu_)

    bo_ = Nucleus.Collection.is_in(n1_, n2_)

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

    ax = um * 0.008

    Nucleus.Plotly.writ(
        ht,
        (
            merge(
                tr,
                Dict(
                    "y" => nu_[i1_],
                    "x" => xc_[i1_],
                    "text" => n1_[i1_],
                    "fillcolor" => Nucleus.Color.BL,
                ),
            ),
            merge(
                tr,
                Dict(
                    "y" => nu_[i2_],
                    "x" => xc_[i2_],
                    "text" => n1_[i2_],
                    "fillcolor" => Nucleus.Color.RE,
                ),
            ),
            Dict(
                "yaxis" => "y2",
                "y" => zeros(sum(bo_)),
                "x" => xc_[bo_],
                "text" => n1_[bo_],
                "mode" => "markers",
                "marker" => Dict(
                    "symbol" => "line-ns",
                    "size" => 24,
                    "line" => Dict("width" => 2, "color" => "#000000cc"),
                ),
                "hoverinfo" => "x+text",
            ),
            merge(
                tr,
                Dict(
                    "yaxis" => "y3",
                    "y" => cu_,
                    "x" => xc_,
                    "text" => n1_,
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
                            "x" => 1 - ax,
                            "xanchor" => "right",
                            "text" => a2,
                            "font" => Dict("color" => Nucleus.Color.RE),
                        ),
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => um + ax,
                            "xanchor" => "left",
                            "text" => a1,
                            "font" => Dict("color" => Nucleus.Color.BL),
                        ),
                    ),
                ),
            ),
            la,
        ),
    )

end

function writ(pr, al, n1_, N, n3_, n2__, n4_, E, um = 2, la = Dict{String, Any}(); ke_...)

    i1_ = findall(en_ -> all(!isnan, en_), eachrow(E))

    n3_ = n3_[i1_]

    n2__ = n2__[i1_]

    E = E[i1_, :]

    Nucleus.Table.writ("$pr.tsv", Nucleus.Table.make("Set", n3_, n4_, E))

    Nucleus.HeatPlot.writ("$pr.html", n3_, n4_, E, la)

    for i2_ in CartesianIndices(E)[Nucleus.Extreme.index(vec(E), um)]

        i1, n2 = Tuple(i2_)

        n3 = n3_[i1]

        writ(
            "$pr.$(Nucleus.Numbe.text_2(E[i2_])).$(n4_[n2]).$n3.html",
            al,
            n1_,
            N[:, n2],
            n2__[i1],
            Dict("title" => Dict("text" => n3));
            ke_...,
        )

    end

end

end
