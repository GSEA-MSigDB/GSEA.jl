using Test: @test

include("_.jl")

# ---- #

for (st, re) in zip(("ks0", "a0", "da2", "da2w", "da2w0w"), AL_)

    @test GSEA.Algorithm.make(st) === re

end
