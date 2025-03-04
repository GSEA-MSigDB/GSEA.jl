using GSEA

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
