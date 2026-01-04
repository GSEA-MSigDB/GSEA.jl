module Normalization

using StatsBase: mean

using Public

using ..GSEA

function make(en, m1, m2)

    en / (en < 0 ? -m1 : m2)

end

function make!(en_, R)

    for i1 in axes(R, 1)

        ne_, po_ = Public.number_sign(R[i1, :])

        m1 = mean(ne_)

        m2 = mean(po_)

        en_[i1] = make(en_[i1], m1, m2)

        for i2 in axes(R, 2)

            R[i1, i2] = make(R[i1, i2], m1, m2)

        end

    end

end

end
