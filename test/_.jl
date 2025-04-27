using Nucleus

using GSEA

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

const S1_ = ["Aa", "22", "33", "44", "55", "66", "77", "88", "99", "Xx", "Jj", "Qq", "Kk"]

const N1_ = [-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6.0]

const S2_ = ["Aa", "Kk", "Jo"]

const DA = pkgdir(GSEA, "data")

const S3_, N2_ = eachcol(Nucleus.Table.rea(joinpath(DA, "myc.tsv"); select = [1, 2]))

const S4_ =
    GSEA.File.read_gmt(joinpath(DA, "c2.all.v7.1.symbols.gmt"))["COLLER_MYC_TARGETS_UP"]

const DI = GSEA.File.read_gmt(joinpath(DA, "h.all.v7.1.symbols.gmt"))

const S5_ = collect(keys(DI))

const ST__ = collect(values(DI))
