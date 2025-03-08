using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

for na in (
    "Algorithm",
    "CommandLineInterface",
    "Enrichment",
    "File",
    "Interface",
    "Normalization",
    "Plot",
    "Rando",
    "Result",
    "Sort",
)

    @info "🎬 Testing $na"

    run(`julia --project $na.jl`)

end
