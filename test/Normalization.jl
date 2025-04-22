using Test: @test

using GSEA

# ---- #

for (en, m1, m2, re) in ((-1, -2, 2, -0.5), (1, 2, 2, 0.5))

    @test GSEA.Normalization.make(en, m1, m2) === re

end

# ---- #

# 275.084 μs (1500 allocations: 2.34 MiB)

for (en_, R) in ((randn(100), randn(100, 1000)),)

    #@btime GSEA.Normalization.make!($en_, $R)

end
