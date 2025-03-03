module Interface

using Nucleus

function update_sort(fe_, sc_)

    id_ = findall(!isnan, sc_)

    fe_ = fe_[id_]

    sc_ = sc_[id_]

    sortperm!(id_, sc_; rev = true)

    fe_[id_], sc_[id_]

end

function make(al, fe_, sc_, me___; ex = 1.0, mi = 1, ma = 1000, fr = 0.0)

    fe_, sc_ = update_sort(fe_, sc_)

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

function make!(di, al, fe_, sc, ne, se_me_, ns, sa_; st = 0.0, up = 2, ke_ar...)

    if !iszero(st)

        foreach(sc_ -> Nucleus.Normalization.standardize_clamp!(sc_, st), eachcol(sc))

    end

    se_ = collect(keys(se_me_))

    me___ = collect(values(se_me_))

    en = stack((make(al, fe_, sc_, me___; ke_ar...) for sc_ in eachcol(sc)))

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
            update_sort(fe_, sc[:, ia])...,
            me___[is];
            ns = sa,
            la = Dict("title" => Dict("text" => se)),
        )

    end

    se_, en

end

end
