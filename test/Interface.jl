using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Nucleus

include("_.jl")

# ---- #

# 3.419 ms (35 allocations: 1.94 MiB)
# 2.975 ms (35 allocations: 1.94 MiB)
# 21.604 ms (35 allocations: 1.94 MiB)
# 21.615 ms (35 allocations: 1.94 MiB)
# 13.392 ms (35 allocations: 1.94 MiB)
# 12.282 ms (35 allocations: 1.94 MiB)

const N3_ = collect(keys(D1))

const N2__ = collect(values(D1))

for al in AL_

    en_ = GSEA.Interface.make(al, GE_, EX_, N2__)

    #@btime GSEA.Interface.make($al, GE_, EX_, N2__)

    @test !issorted(en_)

    @test N3_[sortperm(en_)][(end - 1):end] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

end
