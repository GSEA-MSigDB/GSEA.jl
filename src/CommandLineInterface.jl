module CommandLineInterface

using Comonicon: @cast, @main

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample, std

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
        Nucleus.Table.make(na_[1], an[!, 1], na_[2:end], Matrix(an[!, 2:end])),
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

function make!(N, st)

    foreach(nu_ -> Nucleus.Normalization.update_0_clamp!(nu_, st), eachcol(N))

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
  - `--maximum`: = 1000. The maximum set size.`
  - `--fraction`: = 0. The minimum fraction of set genes present.
  - `--number-of-plots`: = 2.
  - `--more-plots`: = "". ;-separated set names.
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

        make!(N, standard_deviation)

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

function make_normalized(
    ::Union{GSEA.Algorithm.KS0, GSEA.Algorithm.A0},
    en,
    m1,
    m2,
    ::Real,
    ::Real,
)

    en / (en < 0 ? -m1 : m2)

end

function make_normalized(::Any, en, m1, m2, s1, s2)

    if en < 0

        m3 = m1

        s3 = -s1

    else

        m3 = m2

        s3 = s2

    end

    # TODO: Check.
    1 + (en - m3) / (s3 * 3)

end

function make_normalized!(al, en_, R)

    for i1 in axes(R, 1)

        r1_, r2_ = Nucleus.Numbe.ge(R[i1, :])

        m1 = mean(r1_)

        m2 = mean(r2_)

        s1 = std(r1_)

        s2 = std(r2_)

        en_[i1] = make_normalized(al, en_[i1], m1, m2, s1, s2)

        for i2 in axes(R, 2)

            R[i1, i2] = make_normalized(al, R[i1, i2], m1, m2, s1, s2)

        end

    end

end

function make_random(u1, se, al, n1_, nu_, n2__; ke_...)

    R = Matrix{Float64}(undef, lastindex(n2__), u1)

    u2_ = map(lastindex, n2__)

    seed!(se)

    @showprogress for id in 1:u1

        R[:, id] = GSEA.Interface.make(
            al,
            n1_,
            nu_,
            map(um -> sample(n1_, um; replace = false), u2_);
            ke_...,
        )

    end

    R

end

function writ(di, al, n1_, nu_, n3_, n2__, en_, R, um, n4_, a1, a2)

    id_ = findall(!isnan, en_)

    n3_ = n3_[id_]

    n2__ = n2__[id_]

    en_ = en_[id_]

    R = R[id_, :]

    E = Matrix{Float64}(undef, lastindex(id_), 4)

    E[:, 1] = en_

    make_normalized!(al, en_, R)

    E[:, 2] = en_
    @info "" E

    id_, pv_, qv_ = Nucleus.Significance.make(en_, R)

    E[id_, 3] = pv_

    E[id_, 4] = qv_

    Nucleus.Table.writ(
        joinpath(di, "result.tsv"),
        Nucleus.Table.make(
            "Set",
            n3_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
            E,
        ),
    )

    for id in unique!(
        vcat(Nucleus.Extreme.index(E[:, 1], um), filter!(!isnothing, indexin(n4_, n3_))),
    )

        se = n3_[id]

        GSEA.Plot.writ(
            joinpath(di, "$(Nucleus.Numbe.text(en_[id])).$se.html"),
            al,
            n1_,
            nu_,
            n2__[id],
            Dict("title" => Dict("text" => se));
            a1,
            a2,
        )

    end

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
  - `--maximum`: = 1000. The maximum set size.`
  - `--fraction`: = 0. The minimum fraction of set genes present.
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

    al = GSEA.Algorithm.make(algorithm)

    n1_, nu_ = eachcol(Nucleus.Table.rea(tsv; select = [1, 2]))

    di = Nucleus.Dictionary.rea(json)

    n2__ = collect(values(di))

    ke_ = (mi = minimum, ma = maximum, fr = fraction)

    writ(
        directory,
        al,
        n1_,
        nu_,
        collect(keys(di)),
        n2__,
        GSEA.Interface.make(al, n1_, nu_, n2__; ke_...),
        make_random(number_of_permutations, seed, al, n1_, nu_, n2__; ke_...),
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
  - `--maximum`: = 1000. The maximum set size.`
  - `--fraction`: = 0. The minimum fraction of set genes present.
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

    y1 = Nucleus.Table.rea(tsv1)

    y2 = Nucleus.Table.rea(tsv2)

    n1_ = convert(BitVector, collect(y1[1, 2:end]))

    a1_ = y2[!, 1]

    N = Matrix(y2[!, indexin(names(y1)[2:end], names(y2))])

    if !iszero(standard_deviation)

        make!(N, standard_deviation)

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
        Nucleus.Table.make("Feature", a1_, [metric], reshape(n2_, :, 1)),
    )

    al = GSEA.Algorithm.make(algorithm)

    di = Nucleus.Dictionary.rea(json)

    a2__ = collect(values(di))

    ke_ = (mi = minimum, ma = maximum, fr = fraction)

    if permutation == "set"

        R = make_random(number_of_permutations, seed, al, a1_, n2_, a2__; ke_...)

    elseif permutation == "sample"

        R = Matrix{Float64}(undef, lastindex(a2__), number_of_permutations)

        seed!(seed)

        @showprogress for id in 1:number_of_permutations

            R[:, id] = GSEA.Interface.make(
                al,
                a1_,
                map(n2_ -> Nucleus.PairMap.make(fu, shuffle!(n1_), n2_), eachrow(N)),
                a2__;
                ke_...,
            )

        end

    end

    writ(
        directory,
        al,
        a1_,
        n2_,
        collect(keys(di)),
        a2__,
        GSEA.Interface.make(al, a1_, n2_, a2__; ke_...),
        R,
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
