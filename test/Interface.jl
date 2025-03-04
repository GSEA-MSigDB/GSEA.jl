using Test: @test

using GSEA

# ----------------------------------------------------------------------------------------------- #

using Nucleus

# ---- #

const DI = pkgdir(GSEA, "data")

# ---- #

const TE = joinpath(tempdir(), "GSEA")

rm(TE; recursive = true, force = true)

mkdir(TE)

# ---- #
