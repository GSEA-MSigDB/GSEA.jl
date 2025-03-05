using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

for na in (
    "Algorithm",
    #"CommandLineInterface",
    "File",
    "Interface",
    "Plot",
)

    @info "🎬 Testing $na"

    run(`julia --project $na.jl`)

end
