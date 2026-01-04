module Rando

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: sample

using Public

using ..GSEA

function make(um, se, al, s1_, nu_, st__; ke_...)

    R = Matrix{Float64}(undef, lastindex(st__), um)

    um_ = map(s2_ -> lastindex(intersect(s1_, s2_)), st__)

    seed!(se)

    @showprogress for nd in 1:um

        R[:, nd] = GSEA.Interface.make(
            al,
            s1_,
            nu_,
            map(um -> sample(s1_, um; replace = false), um_);
            ke_...,
        )

    end

    R

end

function make(um, se, al, st_, fu, bo_, N, st__; ke_...)

    R = Matrix{Float64}(undef, lastindex(st__), um)

    seed!(se)

    @showprogress for nd in 1:um

        R[:, nd] = GSEA.Interface.make(
            al,
            st_,
            map(nu_ -> Public.make_2(fu, shuffle!(bo_), nu_), eachrow(N)),
            st__;
            ke_...,
        )

    end

    R

end

end
