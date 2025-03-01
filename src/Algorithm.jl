module Algorithm

using Nucleus

struct KS end

struct KSa end

struct KLioM end

struct KLioP end

struct KLi end

struct KLi1 end

function text(al)

    string(al)[6:(end - 2)]

end

function make_normalizer(::Union{KS, KSa}, nu_, ex, bo_)

    s0 = s1 = 0.0

    for id in eachindex(nu_)

        if bo_[id]

            s1 += Nucleus.Numbe.make_exponential(nu_[id], ex)

        else

            s0 += 1.0

        end

    end

    -inv(s0), inv(s1)

end

function make_normalizer(::Any, nu_, ex, bo_)

    s1 = s2 = 0.0

    for id in eachindex(nu_)

        s2 += ab = Nucleus.Numbe.make_exponential(nu_[id], ex)

        if bo_[id]

            s1 += ab

        end

    end

    inv(s1), inv(s2)

end

function make_normalizer(o1, o2)

    inv(inv(o2) - inv(o1))

end

function make!(al::KS, nu_, ex, bo_, cu_)

    o0, o1 = make_normalizer(al, nu_, ex, bo_)

    c2 = a2 = c1 = 0.0

    for id in eachindex(nu_)

        c1 += bo_[id] ? Nucleus.Numbe.make_exponential(nu_[id], ex) * o1 : o0

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

function make!(al::KSa, nu_, ex, bo_, cu_)

    o0, o1 = make_normalizer(al, nu_, ex, bo_)

    cu = su = 0.0

    for id in eachindex(nu_)

        su += cu += bo_[id] ? Nucleus.Numbe.make_exponential(nu_[id], ex) * o1 : o0

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

# TODO: Clip.
const ON = 1.0 + 1e-13

function make!(al::KLioM, nu_, ex, bo_, cu_)

    o1, o2 = make_normalizer(al, nu_, ex, bo_)

    o0 = make_normalizer(o1, o2)

    r0 = r1 = r2 = eps()

    l0 = l1 = l2 = ON

    p0 = p1 = p2 = su = 0.0

    for id in eachindex(nu_)

        ab = Nucleus.Numbe.make_exponential(nu_[id], ex)

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

function make!(al::KLioP, nu_, ex, bo_, cu_)

    o1, o2 = make_normalizer(al, nu_, ex, bo_)

    o0 = make_normalizer(o1, o2)

    r0 = r1 = r2 = eps()

    l0 = l1 = l2 = ON

    p0 = p1 = p2 = su = 0.0

    for id in eachindex(nu_)

        ab = Nucleus.Numbe.make_exponential(nu_[id], ex)

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

        r0 += p0 = d0

        r1 += p1 = d1

        r2 += p2 = d2

        su +=
            cu =
                Nucleus.Information.make_symmetric_kullback_leibler_divergence(r1, r0, r2) -
                Nucleus.Information.make_symmetric_kullback_leibler_divergence(l1, l0, l2)

        if !isnothing(cu_)

            cu_[id] = cu

        end

    end

    su / lastindex(nu_)

end

function make!(al::KLi, nu_, ex, bo_, cu_)

    o1, o2 = make_normalizer(al, nu_, ex, bo_)

    r1 = r2 = eps()

    l1 = l2 = ON

    p1 = p2 = su = 0.0

    for id in eachindex(nu_)

        ab = Nucleus.Numbe.make_exponential(nu_[id], ex)

        d1 = bo_[id] ? ab * o1 : 0.0

        d2 = ab * o2

        r1 += d1

        r2 += d2

        l1 -= p1

        l2 -= p2

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

        p1 = d1

        p2 = d2

    end

    su / lastindex(nu_)

end

function make!(al::KLi1, nu_, ex, bo_, cu_)

    um = lastindex(nu_)

    o1, _ = make_normalizer(al, nu_, ex, bo_)

    d2 = inv(um)

    r1 = r2 = eps()

    l1 = ON

    l2 = ON + d2

    p1 = su = 0.0

    for id in eachindex(nu_)

        d1 = bo_[id] ? Nucleus.Numbe.make_exponential(nu_[id], ex) * o1 : 0.0

        r1 += d1

        r2 += d2

        l1 -= p1

        l2 -= d2

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

        p1 = d1

    end

    su / um

end

end
