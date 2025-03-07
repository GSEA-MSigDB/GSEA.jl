module Plot

using Printf: @sprintf

using Nucleus

using ..GSEA

function writ(ht, al, t1_, y1_, t2_, la = Dict{String, Any}(); a1 = "Low", a2 = "High")

    t1_, y1_ = GSEA.Sort.make(t1_, y1_)

    um = lastindex(t1_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    xc_ = collect(1:um)

    i1_ = findall(<(0.0), y1_)

    i2_ = findall(>=(0.0), y1_)

    bo_ = Nucleus.Collection.is_in(t1_, t2_)

    y3_ = Vector{Float64}(undef, um)

    en = GSEA.Enrichment.make!(al, y1_, bo_, y3_)

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
                    "y" => y1_[i1_],
                    "x" => xc_[i1_],
                    "text" => t1_[i1_],
                    "fillcolor" => Nucleus.Color.BL,
                ),
            ),
            merge(
                tr,
                Dict(
                    "y" => y1_[i2_],
                    "x" => xc_[i2_],
                    "text" => t1_[i2_],
                    "fillcolor" => Nucleus.Color.RE,
                ),
            ),
            Dict(
                "yaxis" => "y2",
                "y" => zeros(sum(bo_)),
                "x" => xc_[bo_],
                "text" => t1_[bo_],
                "mode" => "markers",
                "marker" => Dict(
                    "symbol" => "line-ns",
                    "size" => 24,
                    "line" => Dict(
                        "width" => 2,
                        "color" => Nucleus.Color.make(Nucleus.Color.S2, 0.72),
                    ),
                ),
                "hoverinfo" => "x+text",
            ),
            merge(
                tr,
                Dict(
                    "yaxis" => "y3",
                    "y" => y3_,
                    "x" => xc_,
                    "text" => t1_,
                    "fillcolor" => "#07fa07",
                ),
            ),
        ),
        Nucleus.Dictionary.make(
            Dict(
                "showlegend" => false,
                "yaxis3" =>
                    Dict("domain" => (0.328, 1), "title" => Dict("text" => "Δ Enrichment")),
                "yaxis2" => Dict(
                    "domain" => (0.248, 0.32),
                    "title" => Dict("text" => "Set"),
                    "tickvals" => (),
                ),
                "yaxis" => Dict("domain" => (0, 0.24), "title" => Dict("text" => "Score")),
                "xaxis" => Dict(
                    "title" => Dict("text" => "Feature ($um)"),
                    "showspikes" => true,
                    "spikemode" => "across",
                    "spikedash" => "solid",
                    "spikethickness" => 0.8,
                    "spikecolor" => "#000000",
                ),
                "annotations" => (
                    Dict(
                        "yref" => "paper",
                        "xref" => "paper",
                        "y" => 1.064,
                        "text" => "Enrichment = <b>$(@sprintf "%.4g" en)</b>",
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

function writ(pr, al, n1_, N, n3_, n2__, n4_, E, um = 2; ke_...)

    i1_ = findall(en_ -> all(!isnan, en_), eachrow(E))

    n3_ = n3_[i1_]

    n2__ = n2__[i1_]

    E = E[i1_, :]

    n4 = "Set"

    Nucleus.Table.writ("$pr.tsv", Nucleus.Table.make(n4, n3_, n4_, E))

    Nucleus.HeatPlot.writ(
        "$pr.html",
        n3_,
        n4_,
        E,
        Dict("yaxis" => Dict("title" => n4), "xaxis" => Dict("title" => "Sample")),
    )

    i2_ = findall(!isnan, E)

    for i3_ in CartesianIndices(E)[i2_][Nucleus.Extreme.index(E[i2_], um),]

        i1, n2 = Tuple(i3_)

        n3 = n3_[i1]

        writ(
            "$pr.$(Nucleus.Numbe.text(E[i3_])).$(n4_[n2]).$n3.html",
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
