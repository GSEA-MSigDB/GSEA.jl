module CommandLineInterface

using Comonicon: @cast, @main

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

using Nucleus

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
