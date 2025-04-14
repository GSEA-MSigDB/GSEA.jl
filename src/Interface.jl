module Interface

using Nucleus

using ..GSEA

function make(al, st_, nu_, st__; u1 = 1, u2 = 1000, pr = 0)

    st_, nu_ = GSEA.Sort.make(st_, nu_)

    di = Dict(st_[id] => id for id in eachindex(st_))

    bo_ = falses(length(st_))

    en_ = Vector{Float64}(undef, lastindex(st__))

    for id in eachindex(st__)

        st_ = st__[id]

        Nucleus.Collection.is_in!(bo_, di, st_)

        um = sum(bo_)

        en_[id] =
            um < u1 || u2 < um || um / lastindex(st_) < pr ? NaN :
            GSEA.Enrichment.make!(al, nu_, bo_)

        fill!(view(bo_, bo_), false)

    end

    en_

end

end
