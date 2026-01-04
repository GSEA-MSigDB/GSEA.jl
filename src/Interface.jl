module Interface

using ..GSEA

function make(al, s1_, nu_, st__; mi = 1, ma = 1000, pr = 0)

    s1_, nu_ = GSEA.Sort.make(s1_, nu_)

    di = Dict(s1_[nd] => nd for nd in eachindex(s1_))

    bo_ = fill(false, length(s1_))

    en_ = Vector{Float64}(undef, lastindex(st__))

    for nd in eachindex(st__)

        s2_ = st__[nd]

        GSEA.is_in!(bo_, di, s2_)

        um = sum(bo_)

        en_[nd] =
            um < mi || ma < um || um / lastindex(s2_) < pr ? NaN :
            GSEA.Enrichment.make!(al, nu_, bo_)

        fill!(bo_, false)

    end

    en_

end

end
