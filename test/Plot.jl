using GSEA

include("_.jl")

# ---- #

for (n1_, nu_, n2_) in
    ((C1_, LI_, C2_), (C1_, LI_, C3_), (GE_, EX_, D2["COLLER_MYC_TARGETS_UP"])),
    al in AL_

    GSEA.Plot.writ("", al, n1_, nu_, n2_, Dict("title" => Dict("text" => string(al))))

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
        Dict{String, Any}(
            "title" => Dict("text" => string(al)),
            "yaxis" => Dict("title" => Dict("text" => 'H')),
        );
        a2 = "High Gene Expression",
    )

end
