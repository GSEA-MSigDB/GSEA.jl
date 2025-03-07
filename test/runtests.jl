using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

for na in
    ("Algorithm", "CommandLineInterface", "Enrichment", "File", "Interface", "Plot", "Sort")

    @info "🎬 Testing $na"

    run(`julia --project $na.jl`)

end
