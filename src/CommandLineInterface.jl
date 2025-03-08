module CommandLineInterface

using Comonicon: @cast, @main

using Nucleus

using ..GSEA

"""
Convert .cls to .tsv.

# Arguments

  - `tsv`:
  - `cls`:
"""
@cast function cls(tsv, cls)

    an = GSEA.File.read_cls(cls)

    na_ = names(an)

    Nucleus.Table.writ(
        tsv,
        Nucleus.Table.make(na_[1], an[1, 1], na_[2:end], Matrix(an[!, 2:end])),
    )

end

"""
Convert .gct to .tsv.

# Arguments

  - `tsv`:
  - `gct`:
"""
@cast function gct(tsv, gct)

    an = GSEA.File.read_gct(gct)

    Nucleus.Table.writ(
        tsv,
        Nucleus.Table.make("Feature", an[!, 1], names(an)[2:end], Matrix(an[!, 2:end])),
    )

end

"""
Merge .gmts into .json.

# Arguments

  - `json`:
  - `gmts`:
"""
@cast function gmt(json, gmts...)

    Nucleus.Dictionary.writ(json, reduce(merge!, (GSEA.File.read_gmt(gm) for gm in gmts)))

end

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `directory`:
  - `tsv`:
  - `json`:

# Options

  - `--standard-deviation`: = 0. For column-wise normalization. 0 skips normalization.
  - `--algorithm`: = "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: = 1. The minimum set size.
  - `--maximum`: = 1000. The maximum set size.
  - `--fraction`: = 0. The minimum fraction of set members present.
  - `--number-of-plots`: = 2.
  - `--low`: = "Low".
  - `--high`: = "High".
"""
@cast function data_rank(
    directory,
    tsv,
    json;
    standard_deviation::Real = 0,
    algorithm = "ks0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Real = 0,
    number_of_plots::Int = 2,
    low = "Low",
    high = "High",
)

    al = GSEA.Algorithm.make(algorithm)

    an = Nucleus.Table.rea(tsv)

    n1_ = an[!, 1]

    N = Matrix(an[!, 2:end])

    if !iszero(standard_deviation)

        foreach(
            nu_ -> Nucleus.Normalization.update_0_clamp!(nu_, standard_deviation),
            eachcol(N),
        )

    end

    di = Nucleus.Dictionary.rea(json)

    n2__ = collect(values(di))

    GSEA.Plot.writ(
        joinpath(directory, "result"),
        al,
        n1_,
        N,
        collect(keys(di)),
        n2__,
        names(an)[2:end],
        hcat(
            (
                GSEA.Interface.make(
                    al,
                    n1_,
                    nu_,
                    n2__;
                    mi = minimum,
                    ma = maximum,
                    fr = fraction,
                ) for nu_ in eachcol(N)
            )...,
        ),
        number_of_plots;
        a1 = low,
        a2 = high,
    )

end

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `directory`:
  - `tsv`:
  - `json`:

# Options

  - `--algorithm`: = "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: = 1. The minimum set size.
  - `--maximum`: = 1000. The maximum set size.
  - `--fraction`: = 0. The minimum fraction of set members present.
  - `--number-of-permutations`: = 100.
  - `--seed`: = 20150603.
  - `--number-of-plots`: = 2.
  - `--more-plots`: = "". ;-separated set names.
  - `--low`: = "Low".
  - `--high`: = "High".
"""
@cast function user_rank(
    directory,
    tsv,
    json;
    algorithm = "ks0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Real = 0,
    number_of_permutations::Int = 100,
    seed::Int = 20150603,
    number_of_plots::Int = 2,
    more_plots = "",
    low = "Low",
    high = "High",
)

    n1_, nu_ = eachcol(Nucleus.Table.rea(tsv; select = [1, 2]))

    al = GSEA.Algorithm.make(algorithm)

    di = Nucleus.Dictionary.rea(json)

    n2__ = collect(values(di))

    ke_ = (mi = minimum, ma = maximum, fr = fraction)

    GSEA.Result.writ(
        directory,
        al,
        n1_,
        nu_,
        collect(keys(di)),
        n2__,
        GSEA.Interface.make(al, n1_, nu_, n2__; ke_...),
        GSEA.Rando.make(number_of_permutations, seed, al, n1_, nu_, n2__; ke_...),
        number_of_plots,
        split(more_plots, ';'),
        low,
        high,
    )

end

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `directory`:
  - `tsv1`:
  - `tsv2`:
  - `json`:

# Options

  - `--standard-deviation`: = 0. For column-wise normalization. 0 skips normalization.
  - `--metric`: = "signal-to-noise-ratio" | "mean-difference" | "log-ratio".
  - `--algorithm`: = "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: = 1. The minimum set size.
  - `--maximum`: = 1000. The maximum set size.
  - `--fraction`: = 0. The minimum fraction of set members present.
  - `--permutation`: = "sample" | "set".
  - `--number-of-permutations`: = 100.
  - `--seed`: = 20150603.
  - `--number-of-plots`: = 2.
  - `--more-plots`: = "". ;-separated set names.
  - `--low`: = "Low".
  - `--high`: = "High".
"""
@cast function metric_rank(
    directory,
    tsv1,
    tsv2,
    json;
    standard_deviation::Real = 0,
    metric = "signal-to-noise-ratio",
    algorithm = "ks0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Real = 0,
    permutation = "sample",
    number_of_permutations::Int = 100,
    seed::Int = 20150603,
    number_of_plots::Int = 2,
    more_plots = "",
    low = "Low",
    high = "High",
)

    a1 = Nucleus.Table.rea(tsv1)

    a2 = Nucleus.Table.rea(tsv2)

    n1_ = convert(BitVector, collect(a1[1, 2:end]))

    m1_ = a2[!, 1]

    N = Matrix(a2[!, indexin(names(a1)[2:end], names(a2))])

    if !iszero(standard_deviation)

        foreach(
            nu_ -> Nucleus.Normalization.update_0_clamp!(nu_, standard_deviation),
            eachcol(N),
        )

    end

    fu = if metric == "mean-difference"

        Nucleus.PairMetric.make_mean_difference

    elseif metric == "log-ratio"

        Nucleus.PairMetric.make_log_ratio

    elseif metric == "signal-to-noise-ratio"

        Nucleus.PairMetric.make_signal_to_noise_ratio

    end

    n2_ = map(n2_ -> Nucleus.PairMap.make(fu, n1_, n2_), eachrow(N))

    Nucleus.Table.writ(
        joinpath(directory, "metric.tsv"),
        Nucleus.Table.make("Feature", m1_, [metric], reshape(n2_, :, 1)),
    )

    al = GSEA.Algorithm.make(algorithm)

    di = Nucleus.Dictionary.rea(json)

    m2__ = collect(values(di))

    ke_ = (mi = minimum, ma = maximum, fr = fraction)

    GSEA.Result.writ(
        directory,
        al,
        m1_,
        n2_,
        collect(keys(di)),
        m2__,
        GSEA.Interface.make(al, m1_, n2_, m2__; ke_...),
        if permutation == "set"

            GSEA.Rando.make(number_of_permutations, seed, al, m1_, n2_, m2__; ke_...)

        elseif permutation == "sample"

            GSEA.Rando.make(
                number_of_permutations,
                seed,
                al,
                m1_,
                fu,
                n1_,
                N,
                m2__;
                ke_...,
            )

        end,
        number_of_plots,
        split(more_plots, ';'),
        low,
        high,
    )

end

"""
"""
@main

end
