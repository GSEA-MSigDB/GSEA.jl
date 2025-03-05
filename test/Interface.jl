using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Nucleus

include("_.jl")

# ---- #

for (na_, nu_, re) in (('a':'f', [1, NaN, 3, NaN, 5], (['e', 'c', 'a'], [5.0, 3.0, 1.0])),)

    @test GSEA.Interface.update(na_, nu_) == re

end

# ---- #

# 3.103 ms (13 allocations: 803.23 KiB)
# 2.656 ms (13 allocations: 803.23 KiB)
# 21.366 ms (13 allocations: 803.23 KiB)
# 21.389 ms (13 allocations: 803.23 KiB)
# 13.123 ms (13 allocations: 803.23 KiB)
# 11.968 ms (13 allocations: 803.23 KiB)

const N3_ = collect(keys(D1))

const N2__ = collect(values(D1))

GE_, EX_ = GSEA.Interface.update(GE_, EX_)

for al in AL_

    en_ = GSEA.Interface.make(al, GE_, EX_, N2__)

    #@btime GSEA.Interface.make($al, GE_, EX_, N2__)

    @test !issorted(en_)

    @test N3_[sortperm(en_)][(end - 1):end] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

end
