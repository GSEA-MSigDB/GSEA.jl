module Rando

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: sample

using Nucleus

using ..GSEA

function make(u1, se, al, n1_, nu_, n2__; ke_...)

    R = Matrix{Float64}(undef, lastindex(n2__), u1)

    u2_ = map(n2_ -> lastindex(intersect(n1_, n2_)), n2__)

    seed!(se)

    @showprogress for id in 1:u1

        R[:, id] = GSEA.Interface.make(
            al,
            n1_,
            nu_,
            map(um -> sample(n1_, um; replace = false), u2_);
            ke_...,
        )

    end

    R

end

function make(u1, se, al, n1_, fu, u1_, N, n2__; ke_...)

    R = Matrix{Float64}(undef, lastindex(n2__), u1)

    seed!(se)

    @showprogress for id in 1:u1

        R[:, id] = GSEA.Interface.make(
            al,
            n1_,
            map(nu_ -> Nucleus.PairMap.make(fu, shuffle!(u1_), nu_), eachrow(N)),
            n2__;
            ke_...,
        )

    end

    R

end

end
