using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

for (na_, nu_, re) in (('a':'f', [1, NaN, 3, NaN, 5], (['e', 'c', 'a'], [5.0, 3.0, 1.0])),)

    @test GSEA.update(na_, nu_) == re

end

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
