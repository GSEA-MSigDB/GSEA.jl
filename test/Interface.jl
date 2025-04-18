using Test: @test

using GSEA

include("_.jl")

# ---- #

# 1.781 ms (35 allocations: 1.97 MiB)
# 1.777 ms (35 allocations: 1.97 MiB)
# 7.483 ms (35 allocations: 1.97 MiB)
# 7.907 ms (35 allocations: 1.97 MiB)
# 12.868 ms (35 allocations: 1.97 MiB)

const ST_ = collect(keys(D1))

const ST__ = collect(values(D1))

for al in AL_

    en_ = GSEA.Interface.make(al, GE_, EX_, ST__)

    #@btime GSEA.Interface.make($al, GE_, EX_, ST__)

    @test ST_[sortperm(en_)][(end - 1):end] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

end
