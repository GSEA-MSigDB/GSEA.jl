module GSEA

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

    include("$na.jl")

end

end
