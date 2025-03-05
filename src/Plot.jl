module Plot

using Printf: @sprintf

using Nucleus

using ..GSEA

function writ(
    ht,
    al,
    n1_,
    nu_,
    n2_,
    la = Dict{String, Any}();
    ex = 1,
    xa = "Feature",
    y1 = "Score",
    a1 = "Low",
    a2 = "High",
)

    um = lastindex(n1_)

    xc_ = collect(1:um)

    bo_ = Nucleus.Collection.is_in(n1_, n2_)

    cu_ = Vector{Float64}(undef, um)

    en = GSEA.Algorithm.make!(al, nu_, ex, bo_, cu_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    i1_ = findall(<(0.0), nu_)

    i2_ = findall(>=(0.0), nu_)

    y2 = "Δ Enrichment"

    tr_ = [
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
                "y" => cu_,
                "x" => xc_,
                "text" => n1_,
                "fillcolor" => "#07fa07",
            ),
        ),
    ]

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
        tr_,
        Nucleus.Dictionary.make(
            Dict(
                "showlegend" => false,
                "yaxis3" => Dict("domain" => (0.328, 1), "title" => Dict("text" => y2)),
                "yaxis2" => Dict(
                    "domain" => (0.248, 0.32),
                    "title" => Dict("text" => "Set"),
                    "tickvals" => (),
                ),
                "yaxis" => Dict("domain" => (0, 0.24), "title" => Dict("text" => y1)),
                "xaxis" => Dict(
                    "title" => Dict("text" => "$xa ($um)"),
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
                        "font" => Dict("size" => 24, "color" => Nucleus.Color.DA),
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

function writ(pr, al, n1_, N, n3_, n2__, E, xc_, ex = 1, um = 2)

    i1_ = findall(en_ -> all(!isnan, en_), eachrow(E))

    n3_ = n3_[i1_]

    n2__ = n2__[i1_]

    E = E[i1_, :]

    Nucleus.Table.writ("$pr.tsv", Nucleus.Table.make("Set", n3_, xc_, E))

    Nucleus.HeatPlot.writ(
        "$pr.html",
        n3_,
        xc_,
        E,
        Dict("yaxis" => Dict("title" => "Set"), "xaxis" => Dict("title" => "Sample")),
    )

    i2_ = findall(!isnan, E)

    for i3_ in CartesianIndices(E)[i2_][Nucleus.Extreme.index(E[i2_], um)]

        i4, i5 = Tuple(i3_)

        n3 = n3_[i4]

        y1 = xc_[i5]

        writ(
            "$pr.$(Nucleus.Numbe.text(E[i3_])).$y1.$n3.html",
            al,
            GSEA.Interface.update(n1_, N[:, i5])...,
            n2__[i4],
            Dict("title" => Dict("text" => n3));
            ex,
            y1,
        )

    end

end

end
