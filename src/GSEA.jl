module GSEA

# ----------------------------------------------------------------------------------------------- #

using Comonicon: @cast, @main

using Printf: @sprintf

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

using Nucleus

function read_cls(cl)

    l1, l2, l3 = readlines(cl)

    l2 = l2[2:end]

    s3_ = split(l3)

    na = "Phenotype"

    na_ = map(id -> "Sample $id", eachindex(s3_))

    if l1 == "#numeric"

        Nucleus.Table.make(na, l2, na_, [parse(Float64, st) for _ in 1:1, st in s3_])

    else

        s1_ = split(l1)

        s2_ = split(l2)

        @assert parse(Int, s1_[1]) == lastindex(s3_)

        @assert parse(Int, s1_[2]) == lastindex(s2_) == lastindex(unique(s3_))

        di = Dict(st => id for (id, st) in enumerate(s2_))

        Nucleus.Table.make(na, join(s2_, '_'), na_, [di[st] for _ in 1:1, st in s3_])

    end

end

function read_gct(gc)

    Nucleus.Table.rea(gc; header = 3, drop = ["Description"])

end

function read_gmt(gm)

    di = Dict{String, Vector{String}}()

    for li in eachline(gm)

        sp_ = split(li, '\t')

        na = sp_[1]

        @assert !haskey(di, na)

        di[na] = filter!(!isempty, sp_[3:end])

    end

    di

end

"""
Convert .cls to .tsv.

# Arguments

  - `tsv`:
  - `cls`:
"""
@cast function cls(tsv, cls)

    an = read_cls(cls)

    na_ = names(an)

    Nucleus.Table.writ(
        tsv,
        Nucleus.Table.make(na_[1], an[!, 1], na_[2:end], Matrix(an[!, 2:end])),
    )

end

"""
Convert .gct to .tsv.

# Arguments

  - `tsv`:
  - `gct`:
"""
@cast function gct(tsv, gct)

    an = read_gct(gct)

    Nucleus.Table.writ(
        tsv,
        Nucleus.Table.make("Feature", an[!, 1], names(an)[2:end], Matrix(an[!, 2:end])),
    )

end

"""
Merge .gmts into .json.

# Arguments

  - `json`:
  - `gmt_`:
"""
@cast function gmt(json, gmt_...)

    Nucleus.Dictionary.writ(json, reduce(merge!, (read_gmt(gm) for gm in gmt_)))

end

struct KS end

struct KSa end

struct KLioM end

struct KLioP end

struct KLi end

struct KLi1 end

function text(al)

    string(al)[6:(end - 2)]

end

function make_normalizer(::Union{KS, KSa}, nu_, ex, bo_)

    s0 = s1 = 0.0

    for id in eachindex(nu_)

        if bo_[id]

            s1 += Nucleus.Numbe.make_exponential(nu_[id], ex)

        else

            s0 += 1.0

        end

    end

    -inv(s0), inv(s1)

end

function make_normalizer(::Any, nu_, ex, bo_)

    s1 = s2 = 0.0

    for id in eachindex(nu_)

        s2 += ab = Nucleus.Numbe.make_exponential(nu_[id], ex)

        if bo_[id]

            s1 += ab

        end

    end

    inv(s1), inv(s2)

end

function make_normalizer(o1, o2)

    inv(inv(o2) - inv(o1))

end

function _enrich!(al::KS, nu_, ex, bo_, cu_)

    o0, o1 = make_normalizer(al, nu_, ex, bo_)

    c2 = a2 = c1 = 0.0

    for id in eachindex(nu_)

        c1 += bo_[id] ? Nucleus.Numbe.make_exponential(nu_[id], ex) * o1 : o0

        if !isnothing(cu_)

            cu_[id] = c1

        end

        a1 = abs(c1)

        if a2 < a1

            a2 = a1

            c2 = c1

        end

    end

    c2

end

function _enrich!(al::KSa, nu_, ex, bo_, cu_)

    o0, o1 = make_normalizer(al, nu_, ex, bo_)

    cu = su = 0.0

    for id in eachindex(nu_)

        su += cu += bo_[id] ? Nucleus.Numbe.make_exponential(nu_[id], ex) * o1 : o0

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

# TODO: Clip.
const ON = 1.0 + 1e-13

function _enrich!(al::KLioM, nu_, ex, bo_, cu_)

    o1, o2 = make_normalizer(al, nu_, ex, bo_)

    o0 = make_normalizer(o1, o2)

    r0 = r1 = r2 = eps()

    l0 = l1 = l2 = ON

    p0 = p1 = p2 = su = 0.0

    for id in eachindex(nu_)

        ab = Nucleus.Numbe.make_exponential(nu_[id], ex)

        if bo_[id]

            d0 = 0.0

            d1 = ab * o1

        else

            d0 = ab * o0

            d1 = 0.0

        end

        d2 = ab * o2

        l0 -= p0

        l1 -= p1

        l2 -= p2

        r0 += p0 = d0

        r1 += p1 = d1

        r2 += p2 = d2

        su +=
            cu =
                Nucleus.Information.make_antisymmetric_kullback_leibler_divergence(
                    r1,
                    r0,
                    r2,
                ) - Nucleus.Information.make_antisymmetric_kullback_leibler_divergence(
                    l1,
                    l0,
                    l2,
                )

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

function _enrich!(al::KLioP, nu_, ex, bo_, cu_)

    o1, o2 = make_normalizer(al, nu_, ex, bo_)

    o0 = make_normalizer(o1, o2)

    r0 = r1 = r2 = eps()

    l0 = l1 = l2 = ON

    p0 = p1 = p2 = su = 0.0

    for id in eachindex(nu_)

        ab = Nucleus.Numbe.make_exponential(nu_[id], ex)

        if bo_[id]

            d0 = 0.0

            d1 = ab * o1

        else

            d0 = ab * o0

            d1 = 0.0

        end

        d2 = ab * o2

        l0 -= p0

        l1 -= p1

        l2 -= p2

        r0 += p0 = d0

        r1 += p1 = d1

        r2 += p2 = d2

        su +=
            cu =
                Nucleus.Information.make_symmetric_kullback_leibler_divergence(r1, r0, r2) -
                Nucleus.Information.make_symmetric_kullback_leibler_divergence(l1, l0, l2)

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

function _enrich!(al::KLi, nu_, ex, bo_, cu_)

    o1, o2 = make_normalizer(al, nu_, ex, bo_)

    r1 = r2 = eps()

    l1 = l2 = ON

    p1 = p2 = su = 0.0

    for id in eachindex(nu_)

        ab = Nucleus.Numbe.make_exponential(nu_[id], ex)

        d1 = bo_[id] ? ab * o1 : 0.0

        d2 = ab * o2

        r1 += d1

        r2 += d2

        l1 -= p1

        l2 -= p2

        su +=
            cu = Nucleus.Information.make_antisymmetric_kullback_leibler_divergence(
                r1,
                l1,
                r2,
                l2,
            )

        if !isnothing(cu_)

            cu_[id] = cu

        end

        p1 = d1

        p2 = d2

    end

    su / lastindex(nu_)

end

function _enrich!(al::KLi1, nu_, ex, bo_, cu_)

    um = lastindex(nu_)

    o1, _ = make_normalizer(al, nu_, ex, bo_)

    d2 = inv(um)

    r1 = r2 = eps()

    l1 = ON

    l2 = ON + d2

    p1 = su = 0.0

    for id in eachindex(nu_)

        d1 = bo_[id] ? Nucleus.Numbe.make_exponential(nu_[id], ex) * o1 : 0.0

        r1 += d1

        r2 += d2

        l1 -= p1

        l2 -= d2

        su +=
            cu = Nucleus.Information.make_antisymmetric_kullback_leibler_divergence(
                r1,
                l1,
                r2,
                l2,
            )

        if !isnothing(cu_)

            cu_[id] = cu

        end

        p1 = d1

    end

    su / um

end

function _get_extreme(cu_)

    mi, ma = extrema(cu_)

    ai = abs(mi)

    aa = abs(ma)

    if isapprox(ai, aa)

        (mi, ma)

    elseif aa < ai

        (mi,)

    else

        (ma,)

    end

end

function plot(
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

        i3_ = map(in(_get_extreme(cu_)), cu_)

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

function _select_sort(fe_, sc_)

    id_ = findall(!isnan, sc_)

    fe_ = fe_[id_]

    sc_ = sc_[id_]

    sortperm!(id_, sc_; rev = true)

    fe_[id_], sc_[id_]

end

function enrich(al, fe_, sc_, me___; ex = 1.0, mi = 1, ma = 1000, fr = 0.0)

    fe_, sc_ = _select_sort(fe_, sc_)

    en_ = Vector{Float64}(undef, lastindex(me___))

    bo_ = falses(lastindex(fe_))

    fe_ie = Dict(fe_[ie] => ie for ie in eachindex(fe_))

    for is in eachindex(me___)

        me_ = me___[is]

        Nucleus.Collection.is_in!(bo_, fe_ie, me_)

        ui = sum(bo_)

        en_[is] =
            ui < mi || ma < ui || ui / lastindex(me_) < fr ? NaN :
            _enrich!(al, sc_, ex, bo_, nothing)

        bo_[bo_] .= false

    end

    en_

end

function make_algorithm(al)

    if al == "ks"

        KS()

    elseif al == "ksa"

        KSa()

    elseif al == "kliom"

        KLioM()

    elseif al == "kliop"

        KLioP()

    elseif al == "kli"

        KLi()

    elseif al == "kli1"

        KLi1()

    end

end

function _separat(se_me_)

    collect(keys(se_me_)), collect(values(se_me_))

end

function data_rank!(di, al, fe_, sc, ne, se_me_, ns, sa_; st = 0.0, up = 2, ke_ar...)

    if !iszero(st)

        foreach(sc_ -> Nucleus.Normalization.standardize_clamp!(sc_, st), eachcol(sc))

    end

    se_, me___ = _separat(se_me_)

    en = stack((enrich(al, fe_, sc_, me___; ke_ar...) for sc_ in eachcol(sc)))

    ig_ = map(en_ -> all(!isnan, en_), eachrow(en))

    se_ = se_[ig_]

    me___ = me___[ig_]

    en = en[ig_, :]

    pr = joinpath(di, "enrichment")

    Nucleus.XSampleFeature.writ(pr, ne, se_, ns, sa_, text(al), en)

    ig_ = map(!isnan, en)

    for id_ in CartesianIndices(en)[ig_][Nucleus.Extreme.ge(en[ig_], up)]

        is, ia = Tuple(id_)

        se = se_[is]

        sa = sa_[ia]

        plot(
            joinpath(di, "$(Nucleus.Numbe.shorten(en[is, ia])).$sa.$se.html"),
            al,
            _select_sort(fe_, sc[:, ia])...,
            me___[is];
            ns = sa,
            la = Dict("title" => Dict("text" => se)),
        )

    end

    se_, en

end

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `output_directory`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--standard-deviation`: = 0.0. For normalization by column. 0.0 skips normalization.
  - `--algorithm`: = "ks". "ks" | "ksa" | "kliom" | "kliop" | "kli" | "kli1".
  - `--exponent`: = 1.0.
  - `--minimum-set-size`: = 1.
  - `--maximum-set-size`: = 1000.
  - `--set-fraction`: = 0.0.
  - `--number-of-sets-to-plot`: = 2.
  - `--set-name`: = "Set".
  - `--sample-name`: = "Sample".
"""
@cast function data_rank(
    output_directory,
    feature_x_sample_x_score_tsv,
    set_features_json;
    standard_deviation::Float64 = 0.0,
    algorithm = "ks",
    exponent::Float64 = 1.0,
    minimum_set_size::Int = 1,
    maximum_set_size::Int = 1000,
    set_fraction::Float64 = 0.0,
    number_of_sets_to_plot::Int = 2,
    set_name = "Set",
    sample_name = "Sample",
)

    ta = Nucleus.Table.rea(feature_x_sample_x_score_tsv)

    data_rank!(
        output_directory,
        make_algorithm(algorithm),
        ta[!, 1],
        Matrix(ta[!, 2:end]),
        set_name,
        Nucleus.Dic.rea(set_features_json),
        sample_name,
        names(ta)[2:end];
        st = standard_deviation,
        up = number_of_sets_to_plot,
        ex = exponent,
        mi = minimum_set_size,
        ma = maximum_set_size,
        fr = set_fraction,
    )

end

function _normalize_enrichment(::Union{KS, KSa}, en, mn, mp, ::Any, ::Any)

    en / (en < 0.0 ? -mn : mp)

end

function _normalize_enrichment(::Any, en, mn, mp, sn, sp)

    if en < 0.0

        me = mn

        st = -sn

    else

        me = mp

        st = sp

    end

    # TODO: Check `3 * st` or `(3 * st)`.
    1.0 + (en - me) / 3.0 * st

end

function _normalize_enrichment!(al, en_, R)

    us = lastindex(en_)

    no_ = Vector{Float64}(undef, us)

    for id in 1:us

        en = en_[id]

        ra_ = R[id, :]

        rn_, rp_ = Nucleus.Significance._separate(ra_)

        mn = mean(rn_)

        mp = mean(rp_)

        sn = std(rn_)

        sp = std(rp_)

        no_[id] = _normalize_enrichment(al, en, mn, mp, sn, sp)

        R[id, :] = map(R -> _normalize_enrichment(al, R, mn, mp, sn, sp), ra_)

    end

    no_

end

function _write_plot(di, al, fe_, sc_, ex, se_, me___, en_, R, up, pl_, nf, ns, nl, nh)

    ig_ = map(!isnan, en_)

    se_ = se_[ig_]

    me___ = me___[ig_]

    en_ = en_[ig_]

    R = R[ig_, :]

    no_ = _normalize_enrichment!(al, en_, R)

    pn_, qn_, pp_, qp_ = Nucleus.Significance.ge(R, no_)

    Nucleus.Table.writ(
        joinpath(di, "result.tsv"),
        Nucleus.Table.make(
            "Set",
            se_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
            stack((en_, no_, vcat(pn_, pp_), vcat(qn_, qp_))),
        ),
    )

    fe_, sc_ = _select_sort(fe_, sc_)

    for is in
        unique!(vcat(Nucleus.Extreme.ge(en_, up), filter!(!isnothing, indexin(pl_, se_))))

        se = se_[is]

        plot(
            joinpath(di, "$(Nucleus.Numbe.shorten(en_[is])).$se.html"),
            al,
            fe_,
            sc_,
            me___[is];
            ex,
            nf,
            ns,
            nl,
            nh,
            la = Dict("title" => Dict("text" => se)),
        )

    end

end

function _permute_set(ur, sd, al, fe_, sc_, me___; ke_ar...)

    R = Matrix{Float64}(undef, lastindex(me___), ur)

    if !iszero(ur)

        um_ = map(lastindex, me___)

        seed!(sd)

        @showprogress for id in 1:ur

            R[:, id] = enrich(
                al,
                fe_,
                sc_,
                map(um -> sample(fe_, um; replace = false), um_);
                ke_ar...,
            )

        end

    end

    R

end

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `output_directory`:
  - `feature_x_metric_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--algorithm`: = "ks". "ks" | "ksa" | "kliom" | "kliop" | "kli" | "kli1".
  - `--exponent`: = 1.0.
  - `--minimum-set-size`: = 1.
  - `--maximum-set-size`: = 1000.
  - `--set-fraction`: = 0.0.
  - `--number-of-permutations`: = 100.
  - `--random-seed`: = 20150603.
  - `--number-of-sets-to-plot`: = 2.
  - `--more-sets-to-plot`: = "". ;-separated set names.
  - `--feature-name`: = "Gene".
  - `--score-name`: = "My Score".
  - `--low-text`: = "Low".
  - `--high-text`: = "High".
"""
@cast function user_rank(
    output_directory,
    feature_x_metric_x_score_tsv,
    set_features_json;
    algorithm = "ks",
    exponent::Float64 = 1.0,
    minimum_set_size::Int = 1,
    maximum_set_size::Int = 1000,
    set_fraction::Float64 = 0.0,
    number_of_permutations::Int = 100,
    random_seed::Int = 20150603,
    number_of_sets_to_plot::Int = 2,
    more_sets_to_plot = "",
    feature_name = "Gene",
    score_name = "My Score",
    low_text = "Low",
    high_text = "High",
)

    al = make_algorithm(algorithm)

    ta = Nucleus.Table.rea(feature_x_metric_x_score_tsv)

    fe_, sc_ = _select_sort(ta[!, 1], ta[!, 2])

    se_, me___ = _separat(Nucleus.Dic.rea(set_features_json))

    ke_ar = (ex = exponent, mi = minimum_set_size, ma = maximum_set_size, fr = set_fraction)

    _write_plot(
        output_directory,
        al,
        fe_,
        sc_,
        exponent,
        se_,
        me___,
        enrich(al, fe_, sc_, me___; ke_ar...),
        _permute_set(number_of_permutations, random_seed, al, fe_, sc_, me___; ke_ar...),
        number_of_sets_to_plot,
        split(more_sets_to_plot, ';'),
        feature_name,
        score_name,
        low_text,
        high_text,
    )

end

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `output_directory`:
  - `target_x_sample_x_number_tsv`:
  - `feature_x_sample_x_score_tsv`:
  - `set_features_json`:

# Options

  - `--standard-deviation`: = 0.0. For normalization by column. 0.0 skips normalization.
  - `--algorithm`: = "ks". "ks" | "ksa" | "kliom" | "kliop" | "kli" | "kli1".
  - `--exponent`: = 1.0.
  - `--metric`: = "signal-to-noise-ratio". "mean-difference" | "log-ratio" | "signal-to-noise-ratio".
  - `--minimum-set-size`: = 1.
  - `--maximum-set-size`: = 1000.
  - `--set-fraction`: = 0.0.
  - `--permutation`: = "sample". "sample" | "set".
  - `--number-of-permutations`: = 100.
  - `--random-seed`: = 20150603.
  - `--number-of-sets-to-plot`: = 2.
  - `--more-sets-to-plot`: = "". ;-separated set names.
  - `--feature-name`: = "Gene".
  - `--score-name`: = "Signal-to-Noise Ratio".
  - `--low-text`: = "Low".
  - `--high-text`: = "High".
"""
@cast function metric_rank(
    output_directory,
    target_x_sample_x_number_tsv,
    feature_x_sample_x_score_tsv,
    set_features_json;
    standard_deviation::Float64 = 0.0,
    algorithm = "ks",
    exponent::Float64 = 1.0,
    metric = "signal-to-noise-ratio",
    minimum_set_size::Int = 1,
    maximum_set_size::Int = 1000,
    set_fraction::Float64 = 0.0,
    permutation = "sample",
    number_of_permutations::Int = 100,
    random_seed::Int = 20150603,
    number_of_sets_to_plot::Int = 2,
    more_sets_to_plot = "",
    feature_name = "Gene",
    score_name = "Signal-to-Noise Ratio",
    low_text = "Low",
    high_text = "High",
)

    tt = Nucleus.Table.rea(target_x_sample_x_number_tsv)

    tf = Nucleus.Table.rea(feature_x_sample_x_score_tsv)

    vt_ = convert(BitVector, collect(tt[1, 2:end]))

    fe_ = tf[!, 1]

    s1 = Matrix(tf[!, indexin(names(tt)[2:end], names(tf))])

    if !iszero(standard_deviation)

        foreach(
            s1_ -> Nucleus.Normalization.standardize_clamp!(s1_, standard_deviation),
            eachcol(s1),
        )

    end

    fu = if metric == "mean-difference"

        Nucleus.Target.get_mean_difference

    elseif metric == "log-ratio"

        Nucleus.Target.get_log_ratio

    elseif metric == "signal-to-noise-ratio"

        Nucleus.Target.get_signal_to_noise_ratio

    end

    s2_ = map(s1_ -> Nucleus.Target.go(fu, vt_, s1_), eachrow(s1))

    Nucleus.Table.writ(
        joinpath(output_directory, "metric.tsv"),
        Nucleus.Table.make("Feature", fe_, [metric], reshape(s2_, :, 1)),
    )

    al = make_algorithm(algorithm)

    se_, me___ = _separat(Nucleus.Dic.rea(set_features_json))

    ke_ar = (ex = exponent, mi = minimum_set_size, ma = maximum_set_size, fr = set_fraction)

    if permutation == "set"

        R = _permute_set(number_of_permutations, random_seed, al, fe_, s2_, me___; ke_ar...)

    elseif permutation == "sample"

        R = Matrix{Float64}(undef, lastindex(se_), number_of_permutations)

        if 0 < number_of_permutations

            seed!(random_seed)

            @showprogress for id in 1:number_of_permutations

                R[:, id] = enrich(
                    al,
                    fe_,
                    map(s1_ -> Nucleus.Target.go(fu, shuffle!(vt_), s1_), eachrow(s1)),
                    me___;
                    ke_ar...,
                )

            end

        end

    end

    _write_plot(
        output_directory,
        al,
        fe_,
        s2_,
        exponent,
        se_,
        me___,
        enrich(al, fe_, s2_, me___; ke_ar...),
        R,
        number_of_sets_to_plot,
        split(more_sets_to_plot, ';'),
        feature_name,
        score_name,
        low_text,
        high_text,
    )

end

"""
"""
@main

end
