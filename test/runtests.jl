using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

for st in (
    "Algorithm",
    "Enrichment",
    "File",
    "GSEA1",
    "GSEA2",
    "Interface",
    "Normalization",
    "Plot",
    "Rando",
    "Result",
    "Sort",
)

    @info "🎬 Testing $st"

    run(`julia --project $st.jl`)

end
