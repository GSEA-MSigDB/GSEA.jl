include("_.jl")

# ---- #

function make(al)

    Dict("title" => Dict("text" => "$al"))

end

# ---- #

for (s1_, nu_, s2_) in
    ((S1_, N1_, S2_), (S1_, N1_, ["66", "77", "88", "Jo"]), (S3_, N2_, S4_)),
    al in AL_

    GSEA.Plot.writ("", al, s1_, nu_, s2_, make(al))

end

# ---- #

const N = hcat(N2_, N2_ * 2)

for al in AL_

    GSEA.Plot.writ(
        joinpath(TE, "$al"),
        al,
        map(id -> "Sa $id", axes(N, 2)),
        S3_,
        N,
        S5_,
        ST__,
        hcat((GSEA.Interface.make(al, S3_, nu_, ST__) for nu_ in eachcol(N))...),
        1,
        make(al);
        t2 = "High Expression",
    )

end
