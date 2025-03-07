using Test: @test

using GSEA

include("_.jl")

# ---- #

# 3.420 ms (35 allocations: 1.94 MiB)
# 2.980 ms (35 allocations: 1.94 MiB)
# 12.289 ms (35 allocations: 1.94 MiB)
# 13.302 ms (35 allocations: 1.94 MiB)
# 21.676 ms (35 allocations: 1.94 MiB)

const N3_ = collect(keys(D1))

const N2__ = collect(values(D1))

for al in AL_

    en_ = GSEA.Interface.make(al, GE_, EX_, N2__)

    #@btime GSEA.Interface.make($al, GE_, EX_, N2__)

    @test N3_[sortperm(en_)][(end - 1):end] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

end
