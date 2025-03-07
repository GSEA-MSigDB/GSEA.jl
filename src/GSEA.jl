module GSEA

# ----------------------------------------------------------------------------------------------- #

for na in
    ("Algorithm", "CommandLineInterface", "Enrichment", "File", "Interface", "Plot", "Sort")

    include("$na.jl")

end

end
