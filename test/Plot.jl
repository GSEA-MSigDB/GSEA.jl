using Test: @test

using GSEA

using Nucleus

include("_.jl")

# ---- #

for (t1_, yc_, t2_) in
    ((C1_, LI_, C2_), (C1_, LI_, C3_), (GE_, EX_, D2["COLLER_MYC_TARGETS_UP"])),
    al in AL_

    GSEA.Plot.writ("", al, t1_, yc_, t2_, Dict("title" => Dict("text" => string(al))))

end

# ---- #

const N = hcat(EX_, EX_ * 2)

const N3_ = collect(keys(D1))

const N2__ = collect(values(D1))

for al in AL_

    GSEA.Plot.writ(
        joinpath(TE, string(al)),
        al,
        GE_,
        N,
        N3_,
        N2__,
        map(id -> "Sample $id", axes(N, 2)),
        hcat((GSEA.Interface.make(al, GE_, nu_, N2__) for nu_ in eachcol(N))...),
        1,
    )

end
