module Enrichment

using Nucleus

using ..GSEA

function make_delta(::Union{GSEA.Algorithm.KS0, GSEA.Algorithm.A0}, nu_, bo_)

    s0 = s1 = 0.0

    for id in eachindex(nu_)

        if bo_[id]

            s1 += abs(nu_[id])

        else

            s0 += 1.0

        end

    end

    -inv(s0), inv(s1)

end

function make_delta(
    ::Union{GSEA.Algorithm.DA2, GSEA.Algorithm.DA2W, GSEA.Algorithm.DA2W0W},
    nu_,
    bo_,
)

    s1 = s2 = 0.0

    for id in eachindex(nu_)

        s2 += ab = abs(nu_[id])

        if bo_[id]

            s1 += ab

        end

    end

    inv(s1), inv(s2)

end

function make_delta(d1, d2)

    inv(inv(d2) - inv(d1))

end

function make_eps(po)

    ep = eps()

    po < ep ? ep : po

end

function make!(al::GSEA.Algorithm.KS0, nu_, bo_, cu_ = nothing)

    d0, d1 = make_delta(al, nu_, bo_)

    c2 = a2 = c1 = 0.0

    for id in eachindex(nu_)

        c1 += bo_[id] ? d1 * abs(nu_[id]) : d0

        if !isnothing(cu_)

            cu_[id] = c1

        end

        a1 = abs(c1)

        if a2 < a1

            a2 = a1

            c2 = c1

        end

    end

    c2

end

function make!(al::GSEA.Algorithm.A0, nu_, bo_, cu_ = nothing)

    d0, d1 = make_delta(al, nu_, bo_)

    cu = su = 0.0

    for id in eachindex(nu_)

        su += cu += bo_[id] ? d1 * abs(nu_[id]) : d0

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

function make!(al::GSEA.Algorithm.DA2, nu_, bo_, cu_ = nothing)

    d1, _ = make_delta(al, nu_, bo_)

    r1 = r2 = eps()

    l1 = l2 = 1.0

    e1 = su = 0.0

    l2 += e2 = inv(lastindex(nu_))

    for id in eachindex(nu_)

        l1 = make_eps(l1 - e1)

        l2 = make_eps(l2 - e2)

        r1 += e1 = bo_[id] ? d1 * abs(nu_[id]) : 0.0

        r2 += e2

        su +=
            cu = Nucleus.Information.make_antisymmetric_kullback_leibler_divergence(
                r1,
                l1,
                r2,
                l2,
            )

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

function make!(al::GSEA.Algorithm.DA2W, nu_, bo_, cu_ = nothing)

    d1, d2 = make_delta(al, nu_, bo_)

    r1 = r2 = eps()

    l1 = l2 = 1.0

    e1 = e2 = su = 0.0

    for id in eachindex(nu_)

        ab = abs(nu_[id])

        l1 = make_eps(l1 - e1)

        l2 = make_eps(l2 - e2)

        r1 += e1 = bo_[id] ? d1 * ab : 0.0

        r2 += e2 = d2 * ab

        su +=
            cu = Nucleus.Information.make_antisymmetric_kullback_leibler_divergence(
                r1,
                l1,
                r2,
                l2,
            )

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

function make!(al::GSEA.Algorithm.DA2W0W, nu_, bo_, cu_ = nothing)

    d1, d2 = make_delta(al, nu_, bo_)

    d0 = make_delta(d1, d2)

    r0 = r1 = r2 = eps()

    l0 = l1 = l2 = 1.0

    e0 = e1 = e2 = su = 0.0

    for id in eachindex(nu_)

        ab = abs(nu_[id])

        l0 = make_eps(l0 - e0)

        l1 = make_eps(l1 - e1)

        l2 = make_eps(l2 - e2)

        r0 += e0 = bo_[id] ? 0.0 : d0 * ab

        r1 += e1 = bo_[id] ? d1 * ab : 0.0

        r2 += e2 = d2 * ab

        su +=
            cu =
                Nucleus.Information.make_antisymmetric_kullback_leibler_divergence(
                    r1,
                    r0,
                    r2,
                ) - Nucleus.Information.make_antisymmetric_kullback_leibler_divergence(
                    l1,
                    l0,
                    l2,
                )

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

end
