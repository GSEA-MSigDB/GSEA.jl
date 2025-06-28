using Test: @test

include("_.jl")

# ---- #

const A1 = AL_[1]

const A3 = AL_[3]

const B1_ = Nucleus.Collection.is_in(S1_, S2_)

const D1 = 1 / 12

const D2 = 1 / 42

const UM = 100000

const R1_ = randn(UM)

const R2_ = rand(push!(falses(9), true), UM)

# ---- #

function test(al, nu_, bo_, r1, r2)

    d1, d2 = GSEA.Enrichment.make_delta(al, nu_, bo_)

    #@btime GSEA.Enrichment.make_delta($al, $nu_, $bo_)

    @test isnothing(r1) || d1 === r1

    @test isnothing(r2) || d2 === r2

end

# ---- #

# 5.708 ns (0 allocations: 0 bytes)
# 48.958 μs (0 allocations: 0 bytes)

for (nu_, bo_, r0, r1) in ((N1_, B1_, -1 / 11, D1), (R1_, R2_, nothing, nothing))

    test(A1, nu_, bo_, r0, r1)

end

# ---- #

# 5.250 ns (0 allocations: 0 bytes)
# 55.958 μs (0 allocations: 0 bytes)

for (nu_, bo_, r1, r2) in ((N1_, B1_, D1, D2), (R1_, R2_, nothing, nothing))

    test(A3, nu_, bo_, r1, r2)

end

# ---- #

for (d1, d2, re) in ((D1, D2, 1 / 30),)

    @test GSEA.Enrichment.make_delta(d1, d2) === re

end

# ---- #

const EP = eps()

const PO = EP * 2

for (nu, re) in ((-1, EP), (0, EP), (PO, PO))

    @test GSEA.Enrichment.make_eps(nu) === re

end

# ---- #

S1_, N1_ = GSEA.Sort.make(S1_, N1_)

S3_, N2_ = GSEA.Sort.make(S3_, N2_)

const B2_ = Nucleus.Collection.is_in(S3_, S4_)

# ---- #

function test(al, nu_, bo_, re)

    @test GSEA.Enrichment.make!(al, nu_, bo_) === re

    #@btime GSEA.Enrichment.make!($al, $nu_, $bo_)

end

# ---- #

# 10.343 ns (0 allocations: 0 bytes)
# 22.208 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((N1_, B1_, -0.5000000000000001), (N2_, B2_, 0.7651927829281453))

    test(A1, nu_, bo_, re)

end

# ---- #

# 10.594 ns (0 allocations: 0 bytes)
# 21.583 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((N1_, B1_, -7.686159401251084e-17), (N2_, B2_, 0.41482514169516305))

    test(AL_[2], nu_, bo_, re)

end

# ---- #

# 84.458 ns (0 allocations: 0 bytes)
# 127.166 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((N1_, B1_, 4.099285014000578e-16), (N2_, B2_, 1.2297916337424049))

    test(A3, nu_, bo_, re)

end

# ---- #

# 86.057 ns (0 allocations: 0 bytes)
# 135.208 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((N1_, B1_, 1.964240735875277e-16), (N2_, B2_, 1.1161382190540838))

    test(AL_[4], nu_, bo_, re)

end

# ---- #

# 154.732 ns (0 allocations: 0 bytes)
# 225.208 μs (0 allocations: 0 bytes)

for (nu_, bo_, re) in ((N1_, B1_, 1.0248212535001446e-16), (N2_, B2_, 1.1181841586127337))

    test(AL_[5], nu_, bo_, re)

end
