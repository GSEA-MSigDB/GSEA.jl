module Normalization

using StatsBase: mean, std

using Nucleus

using ..GSEA

function make(::Union{GSEA.Algorithm.KS0, GSEA.Algorithm.A0}, en, m1, m2, ::Real, ::Real)

    en / (en < 0 ? -m1 : m2)

end

function make(
    ::Union{GSEA.Algorithm.DA2, GSEA.Algorithm.DA2W, GSEA.Algorithm.DA2W0W},
    en,
    m1,
    m2,
    s1,
    s2,
)

    if en < 0

        m3 = m1

        s3 = -s1

    else

        m3 = m2

        s3 = s2

    end

    1 + (en - m3) / (s3 * 3)

end

function make!(al, en_, R)

    for i1 in axes(R, 1)

        ne_, po_ = Nucleus.Numbe.ge(R[i1, :])

        m1 = mean(ne_)

        m2 = mean(po_)

        s1 = std(ne_)

        s2 = std(po_)

        en_[i1] = make(al, en_[i1], m1, m2, s1, s2)

        for i2 in axes(R, 2)

            R[i1, i2] = make(al, R[i1, i2], m1, m2, s1, s2)

        end

    end

end

end
