using GSEA

# ------------------------------------ #

for st in (
    "Algorithm",
    "Enrichment",
    "File",
    "Interface",
    "Normalization",
    "Plot",
    "Rando",
    "Result",
    "Sort",
    "GSEA",
)

    @info "🎬 Testing $st"

    run(`julia --project $st.jl`)

end
