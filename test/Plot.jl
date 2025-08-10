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
# TODO

GSEA.Plot.writ

# ---- #

const N = hcat(N2_, N2_ * 2)

for al in AL_

    GSEA.Plot.writ(
        joinpath(TE, "$al.html"),
        al,
        S3_,
        map(nd -> "Sa $nd", axes(N, 2)),
        N,
        S5_,
        ST__,
        reduce(hcat, GSEA.Interface.make(al, S3_, nu_, ST__) for nu_ in eachcol(N)),
        make(al);
        um = 1,
        t3 = "High Expression",
    )

end
