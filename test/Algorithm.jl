using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using StatsBase: Weights, sample

using Nucleus

# ---- #

const DI = pkgdir(GSEA, "data")

# ---- #

const AL_ = GSEA.Algorithm.KS(),
GSEA.Algorithm.KSa(),
GSEA.Algorithm.KLioM(),
GSEA.Algorithm.KLioP(),
GSEA.Algorithm.KLi(),
GSEA.Algorithm.KLi1()

# ---- #

const N1_ = [-2, -1, -0.5, 0, 0, 0.5, 1, 2, 3.4]

const B1_ = [true, false, true, false, true, true, false, false, true]

const N2_ = randn(100000)

const B2_ = sample([false, true], Weights([0.9, 0.1]), lastindex(N2_))

# ---- #

# 70.355 ns (0 allocations: 0 bytes)
# 70.385 ns (0 allocations: 0 bytes)
# 8.291 ns (0 allocations: 0 bytes)
# 19.790 ns (0 allocations: 0 bytes)
# 86.041 μs (0 allocations: 0 bytes)

const NO = -0.25

for (nu_, ex, bo_, re) in (
    (N1_, 0.1, B1_, (NO, 0.24581982412836917)),
    (N1_, 0.5, B1_, (NO, 0.21402570288861142)),
    (N1_, 1, B1_, (NO, 0.15625)),
    (N1_, 2, B1_, (NO, 0.06226650062266501)),
    (N2_, 1, B2_, nothing),
)

    al = GSEA.Algorithm.KS()

    @test isnothing(re) || GSEA.Algorithm.make_normalizer(al, nu_, ex, bo_) === re

    #@btime GSEA.Algorithm.make_normalizer($al, $nu_, $ex, $bo_)

end

# ---- #

# 116.992 ns (0 allocations: 0 bytes)
# 116.949 ns (0 allocations: 0 bytes)
# 7.375 ns (0 allocations: 0 bytes)
# 25.016 ns (0 allocations: 0 bytes)
# 92.833 μs (0 allocations: 0 bytes)

for (nu_, ex, bo_, re) in (
    (N1_, 0.1, B1_, (0.24581982412836917, 0.14006007078470165)),
    (N1_, 0.5, B1_, (0.21402570288861142, 0.12366213677204271)),
    (N1_, 1, B1_, (0.15625, 0.09615384615384615)),
    (N1_, 2, B1_, (0.06226650062266501, 0.04533091568449683)),
    (N2_, 1, B2_, nothing),
)

    al = GSEA.Algorithm.KLioM()

    @test isnothing(re) || GSEA.Algorithm.make_normalizer(al, nu_, ex, bo_) === re

    #@btime GSEA.Algorithm.make_normalizer($al, $nu_, $ex, $bo_)

end

# ---- #

for (o1, o2, re) in ((1 / 3, 0.5, -1.0),)

    @test GSEA.Algorithm.make_normalizer(o1, o2) === re

end

# ---- #

const EP = eps()

for (nu, re) in ((-EP, EP), (0.0, EP))

    @test GSEA.Algorithm.make_eps(nu) === re

end

# ---- #

# 16.784 ns (0 allocations: 0 bytes)
# 16.324 ns (0 allocations: 0 bytes)
# 283.186 ns (0 allocations: 0 bytes)
# 283.039 ns (0 allocations: 0 bytes)
# 157.326 ns (0 allocations: 0 bytes)
# 155.055 ns (0 allocations: 0 bytes)
# 45.333 μs (0 allocations: 0 bytes)
# 37.583 μs (0 allocations: 0 bytes)
# 411.958 μs (0 allocations: 0 bytes)
# 411.875 μs (0 allocations: 0 bytes)
# 246.583 μs (0 allocations: 0 bytes)
# 224.542 μs (0 allocations: 0 bytes)

const NA_, N3_ =
    eachcol(reverse!(Nucleus.Table.rea(joinpath(DI, "myc.tsv"); select = [1, 2])))

const ME_ =
    GSEA.File.read_gmt(joinpath(DI, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]

for (nu_, bo_, re_) in (
        (
            [6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6.0],
            Nucleus.Collection.is_in(
                ['K', 'Q', 'J', 'X', '9', '8', '7', '6', '5', '4', '3', '2', 'A'],
                ['K', 'A'],
            ),
            (-0.5, 0.0, 0.0, 0.0, 0.0, 0.0),
        ),
        (
            N3_,
            Nucleus.Collection.is_in(NA_, ME_),
            (
                0.7651927829281453,
                0.41482514169516305,
                1.1181841586127337,
                1.1140922794954267,
                1.1161382190540838,
                1.2297916337424049,
            ),
        ),
    ),
    (al, re) in zip(AL_, re_)

    ex = 1

    @test isapprox(GSEA.Algorithm.make!(al, nu_, ex, bo_, nothing), re; atol = 1e-15)

    #@btime GSEA.Algorithm.make!($al, $nu_, $ex, $bo_, nothing)

end
