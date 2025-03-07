using Test: @test

using GSEA

include("_.jl")

# ---- #

for (an_, nu_, r1, r2) in (('a':'f', [1, NaN, 3, NaN, 5], ['e', 'c', 'a'], [5, 3, 1.0]),)

    an_, nu_ = GSEA.Sort.make(an_, nu_)

    @test is_egal(an_, r1)

    @test is_egal(nu_, r2)

end
