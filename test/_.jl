using GSEA

using Nucleus

const TE = joinpath(tempdir(), "GSEA")

rm(TE; recursive = true, force = true)

mkdir(TE)

function is_egal(a1_, a2_)

    eltype(a1_) === eltype(a2_) && a1_ == a2_

end

const AL_ = GSEA.Algorithm.KS0(),
GSEA.Algorithm.A0(),
GSEA.Algorithm.DA2(),
GSEA.Algorithm.DA2W(),
GSEA.Algorithm.DA2W0W()

const C1_ = ['A', '2', '3', '4', '5', '6', '7', '8', '9', 'X', 'J', 'Q', 'K']
#const C1_ = ["Aa", "22", "33", "44", "55", "66", "77", "88", "99", "Xx", "Jj", "Qq", "Kk"]

const LI_ = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6.0]

const C2_ = ['A', 'K']
#const C2_ = ["Aa", "Kk"]

const C3_ = ['6', '7', '8']

const DA = pkgdir(GSEA, "data")

const GE_, EX_ = eachcol(Nucleus.Table.rea(joinpath(DA, "myc.tsv"); select = [1, 2]))

const D1 = GSEA.File.read_gmt(joinpath(DA, "h.all.v7.1.symbols.gmt"))

const D2 = GSEA.File.read_gmt(joinpath(DA, "c2.all.v7.1.symbols.gmt"))
