using Test: @test

include("_.jl")

# ---- #

# 1.716 ms (34 allocations: 1.99 MiB)
# 1.695 ms (34 allocations: 1.99 MiB)
# 7.114 ms (34 allocations: 1.99 MiB)
# 7.801 ms (34 allocations: 1.99 MiB)
# 12.673 ms (34 allocations: 1.99 MiB)

for al in AL_

    en_ = GSEA.Interface.make(al, S3_, N2_, ST__)

    #@btime GSEA.Interface.make($al, S3_, N2_, ST__)

    @test S5_[sortperm(en_)][49:50] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

end
