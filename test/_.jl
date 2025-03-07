using GSEA

using Nucleus

# ---- #

const TE = joinpath(tempdir(), "GSEA")

rm(TE; recursive = true, force = true)

mkdir(TE)

# ---- #

function is_egal(a1_, a2_)

    eltype(a1_) === eltype(a2_) && a1_ == a2_

end

# ---- #

const AL_ = GSEA.Algorithm.KS(),
GSEA.Algorithm.KSA(),
GSEA.Algorithm.KLIOM(),
GSEA.Algorithm.KLIOP(),
GSEA.Algorithm.KLI(),
GSEA.Algorithm.KLI1()

# ---- #

const C1_ = ['A', '2', '3', '4', '5', '6', '7', '8', '9', 'X', 'J', 'Q', 'K']

const LI_ = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6.0]

const C2_ = ['A', 'K']

# ---- #

const DI = pkgdir(GSEA, "data")

const GE_, EX_ = eachcol(Nucleus.Table.rea(joinpath(DI, "myc.tsv"); select = [1, 2]))

# ---- #

const D1 = GSEA.File.read_gmt(joinpath(DI, "h.all.v7.1.symbols.gmt"))

const D2 = GSEA.File.read_gmt(joinpath(DI, "c2.all.v7.1.symbols.gmt"))
