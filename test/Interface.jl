using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Nucleus

include("_.jl")

# ---- #

for (na_, nu_, re) in (('a':'f', [1, NaN, 3, NaN, 5], (['e', 'c', 'a'], [5.0, 3.0, 1.0])),)

    @test GSEA.Interface.make_sort(na_, nu_) == re

end

# ---- #

const N1_, NU_ =
    eachcol(reverse!(Nucleus.Table.rea(joinpath(DI, "myc.tsv"); select = [1, 2])))

const IC = GSEA.File.read_gmt(joinpath(DI, "h.all.v7.1.symbols.gmt"))

const N3_ = collect(keys(IC))

const N2__ = collect(values(IC))

# ---- #

# 3.202 ms (32 allocations: 1.86 MiB)
# 2.743 ms (32 allocations: 1.86 MiB)
# 21.579 ms (32 allocations: 1.86 MiB)
# 21.581 ms (32 allocations: 1.86 MiB)
# 13.273 ms (32 allocations: 1.86 MiB)
# 12.055 ms (32 allocations: 1.86 MiB)

for al in AL_

    en_ = GSEA.Interface.make(al, N1_, NU_, N2__)

    @btime GSEA.Interface.make($al, N1_, NU_, N2__)

    @test !issorted(en_)

    n3_, en_ = GSEA.Interface.make_sort(N3_, en_)

    @test n3_[1:2] == ["HALLMARK_MYC_TARGETS_V2", "HALLMARK_MYC_TARGETS_V1"]

end
