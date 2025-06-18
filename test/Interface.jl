using Test: @test

include("_.jl")

# ---- #

# 1.741 ms (34 allocations: 1.99 MiB)
# 1.693 ms (34 allocations: 1.99 MiB)
# 7.299 ms (34 allocations: 1.99 MiB)
# 7.777 ms (34 allocations: 1.99 MiB)
# 12.726 ms (34 allocations: 1.99 MiB)

for al in AL_

    en_ = GSEA.Interface.make(al, S3_, N2_, ST__)

    #@btime GSEA.Interface.make($al, S3_, N2_, ST__)

    @test S5_[sortperm(en_)][49:50] ==
          ["HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2"]

end
