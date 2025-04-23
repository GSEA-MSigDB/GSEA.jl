using GSEA

include("_.jl")

# ---- #

function make(al)

    Dict("title" => Dict("text" => "$al"))

end

# ---- #

for (s1_, nu_, s2_) in
    ((C1_, IN_, C2_), (C1_, IN_, ["66", "77", "88", "Jo"]), (G1_, EX_, G2_)),
    al in AL_

    GSEA.Plot.writ("", al, s1_, nu_, s2_, make(al))

end

# ---- #

const E = hcat(EX_, EX_ * 2)

for al in AL_

    GSEA.Plot.writ(
        joinpath(TE, "$al"),
        al,
        map(id -> "Sa $id", axes(E, 2)),
        G1_,
        E,
        HA_,
        HA__,
        hcat((GSEA.Interface.make(al, G1_, ex_, HA__) for ex_ in eachcol(E))...),
        1,
        make(al);
        t2 = "High Expression",
    )

end
