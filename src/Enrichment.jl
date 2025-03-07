module Enrichment

using Nucleus

using ..GSEA

function make_normalizer(::Union{GSEA.Algorithm.KS0, GSEA.Algorithm.A0}, nu_, bo_)

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

function make_normalizer(::Any, nu_, bo_)

    s1 = s2 = 0.0

    for id in eachindex(nu_)

        s2 += ab = abs(nu_[id])

        if bo_[id]

            s1 += ab

        end

    end

    inv(s1), inv(s2)

end

function make_normalizer(o1, o2)

    inv(inv(o2) - inv(o1))

end

function make_eps(nu)

    ep = eps()

    nu < ep ? ep : nu

end

function make!(al::GSEA.Algorithm.KS0, nu_, bo_, cu_)

    o0, o1 = make_normalizer(al, nu_, bo_)

    c2 = a2 = c1 = 0.0

    for id in eachindex(nu_)

        c1 += bo_[id] ? abs(nu_[id]) * o1 : o0

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

function make!(al::GSEA.Algorithm.A0, nu_, bo_, cu_)

    o0, o1 = make_normalizer(al, nu_, bo_)

    cu = su = 0.0

    for id in eachindex(nu_)

        su += cu += bo_[id] ? abs(nu_[id]) * o1 : o0

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

function make!(al::GSEA.Algorithm.DA2, nu_, bo_, cu_)

    o1, _ = make_normalizer(al, nu_, bo_)

    r1 = r2 = eps()

    l1 = l2 = 1.0

    p1 = su = 0.0

    l2 += d2 = inv(lastindex(nu_))

    for id in eachindex(nu_)

        d1 = bo_[id] ? abs(nu_[id]) * o1 : 0.0

        l1 -= p1

        l2 -= d2

        l1 = make_eps(l1)

        l2 = make_eps(l2)

        r1 += p1 = d1

        r2 += d2

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

function make!(al::GSEA.Algorithm.DA2W, nu_, bo_, cu_)

    o1, o2 = make_normalizer(al, nu_, bo_)

    r1 = r2 = eps()

    l1 = l2 = 1.0

    p1 = p2 = su = 0.0

    for id in eachindex(nu_)

        ab = abs(nu_[id])

        d1 = bo_[id] ? ab * o1 : 0.0

        d2 = ab * o2

        l1 -= p1

        l2 -= p2

        l1 = make_eps(l1)

        l2 = make_eps(l2)

        r1 += p1 = d1

        r2 += p2 = d2

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

function make!(al::GSEA.Algorithm.DA2W0W, nu_, bo_, cu_)

    o1, o2 = make_normalizer(al, nu_, bo_)

    o0 = make_normalizer(o1, o2)

    r0 = r1 = r2 = eps()

    l0 = l1 = l2 = 1.0

    p0 = p1 = p2 = su = 0.0

    for id in eachindex(nu_)

        ab = abs(nu_[id])

        if bo_[id]

            d0 = 0.0

            d1 = ab * o1

        else

            d0 = ab * o0

            d1 = 0.0

        end

        d2 = ab * o2

        l0 -= p0

        l1 -= p1

        l2 -= p2

        l0 = make_eps(l0)

        l1 = make_eps(l1)

        l2 = make_eps(l2)

        r0 += p0 = d0

        r1 += p1 = d1

        r2 += p2 = d2

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
