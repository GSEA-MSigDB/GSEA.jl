using Test: @test

using GSEA

include("_.jl")

# ---- #

for (al, re) in zip(("ks0", "a0", "da2", "da2w", "da2w0w"), AL_)

    @test GSEA.Algorithm.make(al) === re

end
