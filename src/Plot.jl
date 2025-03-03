module Plot

using Comonicon: @cast, @main

using Printf: @sprintf

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

using Nucleus

function writ(
    ht,
    al,
    na_,
    nu_,
    me_,
    la = Dict{String, Any}();
    ex = 1,
    xa = "Feature",
    y1 = "Score",
    a1 = "Low",
    a2 = "High",
)

    um = lastindex(na_)

    xc_ = collect(1:um)

    bo_ = Nucleus.Collection.is_in(na_, me_)

    cu_ = Vector{Float64}(undef, um)

    en = _enrich!(al, nu_, ex, bo_, cu_)

    tr = Dict("mode" => "lines", "line" => Dict("width" => 0), "fill" => "tozeroy")

    i1_ = map(<(0.0), nu_)

    i2_ = map(>=(0.0), nu_)

    y2 = "Δ Enrichment"

    tr_ = [
        merge(
            tr,
            Dict(
                "name" => "- Score",
                "y" => nu_[i1_],
                "x" => xc_[i1_],
                "text" => na_[i1_],
                "fillcolor" => Nucleus.Color.BL,
            ),
        ),
        merge(
            tr,
            Dict(
                "name" => "+ Score",
                "y" => nu_[i2_],
                "x" => xc_[i2_],
                "text" => na_[i2_],
                "fillcolor" => Nucleus.Color.RE,
            ),
        ),
        Dict(
            "yaxis" => "y2",
            "name" => "Set",
            "y" => zeros(sum(bo_)),
            "x" => xc_[bo_],
            "text" => na_[bo_],
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
                "name" => y2,
                "y" => cu_,
                "x" => xc_,
                "text" => na_,
                "fillcolor" => "#07fa07",
            ),
        ),
    ]

    if typeof(al) == KS

        i3_ = map(in(get_extreme(cu_)), cu_)

        push!(
            tr_,
            Dict(
                "yaxis" => "y3",
                "name" => "Extrema",
                "y" => cu_[i3_],
                "x" => xc_[i3_],
                "text" => na_[i3_],
                "mode" => "markers",
                "marker" => Dict(
                    "size" => 32,
                    "color" => Nucleus.Color.make(Nucleus.Color.HU, 0.72),
                ),
            ),
        )

    end

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

end
