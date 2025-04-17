using GSEA

include("_.jl")

# ---- #

function make(al)

    Dict("title" => Dict("text" => "$al"))

end

# ---- #

for (s1_, nu_, s2_) in (
        (C1_, IN_, C2_),
        (C1_, IN_, ["66", "77", "88", "Jo"]),
        (GE_, EX_, D2["COLLER_MYC_TARGETS_UP"]),
    ),
    al in AL_

    GSEA.Plot.writ("", al, s1_, nu_, s2_, make(al))

end

# ---- #

const N = hcat(EX_, EX_ * 2)

const S3_ = collect(keys(D1))

const ST__ = collect(values(D1))

for al in AL_

    GSEA.Plot.writ(
        joinpath(TE, "$al"),
        al,
        map(id -> "Sa $id", axes(N, 2)),
        GE_,
        N,
        S3_,
        ST__,
        hcat((GSEA.Interface.make(al, GE_, nu_, ST__) for nu_ in eachcol(N))...),
        1,
        make(al);
        t2 = "High Gene Expression",
    )

end
