module GSEA

# ----------------------------------------------------------------------------------------------- #

for na in (
    "Algorithm",
    #"CommandLineInterface",
    "File",
    "Interface",
    #"Plot",
)

    include("$na.jl")

end

end
