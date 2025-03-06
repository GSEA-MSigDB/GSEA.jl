using GSEA

const C1_ = ['A', '2', '3', '4', '5', '6', '7', '8', '9', 'X', 'J', 'Q', 'K']

const LI_ = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6.0]

const C2_ = ['A', 'K']

const DI = pkgdir(GSEA, "data")

const GE_, EX_ = eachcol(Nucleus.Table.rea(joinpath(DI, "myc.tsv"); select = [1, 2]))

const D1 = GSEA.File.read_gmt(joinpath(DI, "h.all.v7.1.symbols.gmt"))

const D2 = GSEA.File.read_gmt(joinpath(DI, "c2.all.v7.1.symbols.gmt"))

const AL_ = GSEA.Algorithm.KS(),
GSEA.Algorithm.KSa(),
GSEA.Algorithm.KLioM(),
GSEA.Algorithm.KLioP(),
GSEA.Algorithm.KLi(),
GSEA.Algorithm.KLi1()

const TE = joinpath(tempdir(), "GSEA")

rm(TE; recursive = true, force = true)

mkdir(TE)
