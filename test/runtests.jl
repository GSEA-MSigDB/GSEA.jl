using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

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
    "GSEA1",
    "GSEA2",
)

    @info "🎬 Testing $st"

    run(`julia --project $st.jl`)

end
