module Interface

using Nucleus

using ..GSEA

function make_sort(na_, nu_)

    id_ = findall(!isnan, nu_)

    na_ = na_[id_]

    nu_ = nu_[id_]

    sortperm!(id_, nu_; rev = true)

    na_[id_], nu_[id_]

end

function make(al, n1_, nu_, n2__; ex = 1, mi = 1, ma = 1000, fr = 0)

    n1_, nu_ = make_sort(n1_, nu_)

    en_ = Vector{Float64}(undef, lastindex(n2__))

    bo_ = falses(lastindex(n1_))

    di = Dict(n1_[id] => id for id in eachindex(n1_))

    for id in eachindex(n2__)

        n2_ = n2__[id]

        Nucleus.Collection.is_in!(bo_, di, n2_)

        um = sum(bo_)

        en_[id] =
            um < mi || ma < um || um / lastindex(n2_) < fr ? NaN :
            GSEA.Algorithm.make!(al, nu_, ex, bo_, nothing)

        bo_[bo_] .= false

    end

    en_

end

# TODO: Pick up.
function make!(di, al, n1_, sc, ne, ic, ns, sa_; st = 0.0, up = 2, ke_ar...)

    if !iszero(st)

        foreach(nu_ -> Nucleus.Normalization.standardize_clamp!(nu_, st), eachcol(sc))

    end

    se_ = collect(keys(ic))

    n2__ = collect(values(ic))

    en = stack((make(al, n1_, nu_, n2__; ke_ar...) for nu_ in eachcol(sc)))

    ig_ = findall(en_ -> all(!isnan, en_), eachrow(en))

    se_ = se_[ig_]

    n2__ = n2__[ig_]

    en = en[ig_, :]

    pr = joinpath(di, "enrichment")

    Nucleus.XSampleFeature.writ(pr, ne, se_, ns, sa_, text(al), en)

    ig_ = findall(!isnan, en)

    for id_ in CartesianIndices(en)[ig_][Nucleus.Extreme.ge(en[ig_], up)]

        is, ia = Tuple(id_)

        se = se_[is]

        sa = sa_[ia]

        plot(
            joinpath(di, "$(Nucleus.Numbe.shorten(en[is, ia])).$sa.$se.html"),
            al,
            make_sort(n1_, sc[:, ia])...,
            n2__[is];
            ns = sa,
            la = Dict("title" => Dict("text" => se)),
        )

    end

    se_, en

end

end
