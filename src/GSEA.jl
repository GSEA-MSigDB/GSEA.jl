module GSEA

# ----------------------------------------------------------------------------------------------- #

for na in (
    "Algorithm.jl",
    "CommandLineInterface.jl",
    "File.jl",
    "GSEA.jl",
    "Interface.jl",
    "Plot.jl",
)

    include("$na.jl")

end

end
