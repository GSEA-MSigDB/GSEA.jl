module GSEA

# ----------------------------------------------------------------------------------------------- #

function update(na_, nu_)

    id_ = findall(!isnan, nu_)

    na_ = na_[id_]

    nu_ = nu_[id_]

    sortperm!(id_, nu_; rev = true)

    na_[id_], nu_[id_]

end

for na in (
    "Algorithm",
    #"CommandLineInterface",
    "File",
    "Interface",
    "Plot",
)

    include("$na.jl")

end

end
