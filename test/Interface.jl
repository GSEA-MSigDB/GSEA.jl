using Test: @test

include("_.jl")

# ---- #

# 1.765 ms (35 allocations: 1.97 MiB)
# 1.767 ms (35 allocations: 1.97 MiB)
# 7.232 ms (35 allocations: 1.97 MiB)
# 7.690 ms (35 allocations: 1.97 MiB)
# 12.543 ms (35 allocations: 1.97 MiB)

for al in AL_

    en_ = GSEA.Interface.make(al, S3_, N2_, ST__)

    #@btime GSEA.Interface.make($al, S3_, N2_, ST__)

    @test S5_[sortperm(en_)][49:50] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

end
