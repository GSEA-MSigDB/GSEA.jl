using GSEA

# ------------------------------------ #

for st in (
    "Algorithm",
    "Enrichment",
    "File",
    "Information",
    "Interface",
    "Normalization",
    "Plot",
    "Sort",
    "GSEA",
)

    @info "🎬 Testing $st"

    run(`julia --project $st.jl`)

end
