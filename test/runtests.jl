using GSEA

# ------------------------------------ #

for nd in 1:2

    @info "🎬 Including $nd"

    include("$nd.jl")

end
