module Interface

using Nucleus

using ..GSEA

function update(na_, nu_)

    id_ = findall(!isnan, nu_)

    na_ = na_[id_]

    nu_ = nu_[id_]

    sortperm!(id_, nu_; rev = true)

    na_[id_], nu_[id_]

end

function make(al, n1_, nu_, n2__; ex = 1, mi = 1, ma = 1000, fr = 0)

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

end
