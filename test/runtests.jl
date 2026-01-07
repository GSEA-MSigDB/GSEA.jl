using GSEA

# ------------------------------------ #

for nd in 1:2

    @info "🎬 $nd"

    include("$nd.jl")

end
