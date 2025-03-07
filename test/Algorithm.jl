using Test: @test

using GSEA

include("_.jl")

# ---- #

for (al, re) in zip(("ks", "ksa", "kliom", "kliop", "kli", "kli1"), AL_)

    @test GSEA.Algorithm.make(al) === re

end
