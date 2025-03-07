module Interface

using Nucleus

using ..GSEA

function make(al, n1_, nu_, n2__; mi = 1, ma = 1000, fr = 0)

    n1_, nu_ = GSEA.Sort.make(n1_, nu_)

    bo_ = falses(lastindex(n1_))

    di = Dict(n1_[id] => id for id in eachindex(n1_))

    en_ = Vector{Float64}(undef, lastindex(n2__))

    for id in eachindex(n2__)

        n2_ = n2__[id]

        Nucleus.Collection.is_in!(bo_, di, n2_)

        um = sum(bo_)

        en_[id] =
            um < mi || ma < um || um / lastindex(n2_) < fr ? NaN :
            GSEA.Enrichment.make!(al, nu_, bo_, nothing)

        bo_[bo_] .= false

    end

    en_

end

end
