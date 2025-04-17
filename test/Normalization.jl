using Test: @test

using GSEA

include("_.jl")

# ---- #

for (en, m1, m2, s1, s2, re) in ((-1, -2, 2, 3, 3, -0.5), (1, 2, 2, 3, 3, 0.5))

    @test GSEA.Normalization.make(AL_[1], en, m1, m2, s1, s2) === re

end

# ---- #

for (en, m1, m2, s1, s2, re) in ((-1, -2, 2, 3, 3, nothing), (1, 2, 2, 3, 3, nothing))

    @warn "TODO" GSEA.Normalization.make(AL_[3], en, m1, m2, s1, s2)

end

# ---- #

# 290.208 μs (1500 allocations: 2.34 MiB)
# 247.042 μs (1200 allocations: 2.99 MiB)

for al in (AL_[1], AL_[3])

    #@btime GSEA.Normalization.make!($al, $(randn(100)), $(randn(100, 1000)))

end
