using GSEA

# ------------------------------------ #

for nd in 1:4

    @info "🎬 Testing $nd"

    run(`julia --project $nd.jl`)

end
