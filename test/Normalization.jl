using Test: @test

using GSEA

include("_.jl")

# ---- #

# 383.667 μs (1500 allocations: 2.26 MiB)
# 420.459 μs (1200 allocations: 2.91 MiB)

for al in (AL_[1], AL_[end])

    en_ = randn(100)

    R = randn(100, 1000)

    GSEA.Normalization.make!(al, en_, R)

    # TODO: Test.

    #@btime GSEA.Normalization.make!($al, $en_, $R)

end

# ---- #
# TODO

for () in ()

    GSEA.Normalization.make!

end
