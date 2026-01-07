using GSEA

# ------------------------------------ #

for nd in 1:2

    @info "🎬 $nd"

    run(`julia --project $nd.jl`)

end
