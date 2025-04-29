module Interface

using Nucleus

using ..GSEA

function make(al, s1_, nu_, st__; mi = 1, ma = 1000, pr = 0)

    s1_, nu_ = GSEA.Sort.make(s1_, nu_)

    di = Dict(s1_[id] => id for id in eachindex(s1_))

    bo_ = falses(length(s1_))

    en_ = Vector{Float64}(undef, lastindex(st__))

    for id in eachindex(st__)

        s2_ = st__[id]

        Nucleus.Collection.is_in!(bo_, di, s2_)

        um = sum(bo_)

        en_[id] =
            um < mi || ma < um || um / lastindex(s2_) < pr ? NaN :
            GSEA.Enrichment.make!(al, nu_, bo_)

        fill!(bo_, false)

    end

    en_

end

end
