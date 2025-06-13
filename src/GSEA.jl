module GSEA

# ----------------------------------------------------------------------------------------------- #

for st in (
    "Algorithm",
    "Enrichment",
    "File",
    "Interface",
    "Normalization",
    "Plot",
    "Rando",
    "Result",
    "Sort",
)

    include("$st.jl")

end

using Comonicon: @cast, @main

using Nucleus

"""
Convert .cls to .tsv.

# Arguments

  - `tsv`:
  - `cls`:
"""
@cast function cls(tsv, cls)

    s1, s2, st_, N = File.read_cls(cls)

    Nucleus.Table.writ(tsv, Nucleus.Table.make(s1, s2, st_, N))

end

"""
Convert .gct to .tsv.

# Arguments

  - `tsv`:
  - `gct`:
"""
@cast function gct(tsv, gct)

    st, s1_, s2_, N = File.read_gct(gct)

    Nucleus.Table.writ(tsv, Nucleus.Table.make(st, s1_, s2_, N))

end

"""
Merge .gmts into .json.

# Arguments

  - `json`:
  - `gmts`:
"""
@cast function gmt(json, gmts...)

    Nucleus.Dictionary.writ(json, reduce(merge!, (File.read_gmt(gm) for gm in gmts)))

end

function update!(N, st)

    if iszero(st)

        return

    end

    for nu_ in eachcol(N)

        Nucleus.Normalization.update_0_clamp!(nu_, st)

    end

end

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `directory`:
  - `tsv`:
  - `json`:

# Options

  - `--standard-deviation`: For column-wise normalization. 0 skips normalization.
  - `--exponent`:
  - `--algorithm`: "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
  - `--number-of-plots`:
  - `--low`:
  - `--high`:
"""
@cast function data_rank(
    directory,
    tsv,
    json;
    standard_deviation::Float64 = 0.0,
    exponent::Float64 = 1.0,
    algorithm = "ks0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Float64 = 0.0,
    number_of_plots::Int = 2,
    low = "Low",
    high = "High",
)

    al = Algorithm.make(algorithm)

    st, s1_, s2_, N = Nucleus.Table.ge(Nucleus.Table.rea(tsv))

    update!(N, standard_deviation)

    di = Nucleus.Dictionary.rea(json)

    st__ = collect(values(di))

    Plot.writ(
        joinpath(directory, "result"),
        al,
        s1_,
        s2_,
        N,
        collect(keys(di)),
        st__,
        reduce(
            hcat,
            Interface.make(al, s1_, nu_, st__; mi = minimum, ma = maximum, pr = fraction)
            for nu_ in eachcol(N)
        );
        um = number_of_plots,
        t2 = low,
        t3 = high,
    )

end

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `directory`:
  - `tsv`:
  - `json`:

# Options

  - `--exponent`:
  - `--algorithm`: "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
  - `--number-of-permutations`:
  - `--seed`:
  - `--number-of-plots`:
  - `--more-plots`: ;-separated set names.
  - `--low`:
  - `--high`:
"""
@cast function user_rank(
    directory,
    tsv,
    json;
    exponent::Float64 = 1.0,
    algorithm = "ks0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Float64 = 0.0,
    number_of_permutations::Int = 100,
    seed::Int = 20150603,
    number_of_plots::Int = 2,
    more_plots = "",
    low = "Low",
    high = "High",
)

    st_, nu_ = eachcol(Nucleus.Table.rea(tsv)[!, 1:2])

    al = Algorithm.make(algorithm)

    di = Nucleus.Dictionary.rea(json)

    st__ = collect(values(di))

    ke_ = (mi = minimum, ma = maximum, pr = fraction)

    Result.writ(
        directory,
        al,
        st_,
        nu_,
        collect(keys(di)),
        st__,
        Interface.make(al, st_, nu_, st__; ke_...),
        Rando.make(number_of_permutations, seed, al, st_, nu_, st__; ke_...),
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

  - `--standard-deviation`: For column-wise normalization. 0 skips normalization.
  - `--metric`: "signal-to-noise-ratio" | "mean-difference" | "log-ratio".
  - `--exponent`:
  - `--algorithm`: "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
  - `--permutation`: "sample" | "set".
  - `--number-of-permutations`:
  - `--seed`:
  - `--number-of-plots`:
  - `--more-plots`: ;-separated set names.
  - `--low`:
  - `--high`:
"""
@cast function metric_rank(
    directory,
    tsv1,
    tsv2,
    json;
    standard_deviation::Float64 = 0.0,
    metric = "signal-to-noise-ratio",
    exponent::Float64 = 1.0,
    algorithm = "ks0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Float64 = 0.0,
    permutation = "sample",
    number_of_permutations::Int = 100,
    seed::Int = 20150603,
    number_of_plots::Int = 2,
    more_plots = "",
    low = "Low",
    high = "High",
)

    _, _, s1_, N1 = Nucleus.Table.ge(Nucleus.Table.rea(tsv1))

    st, s2_, s3_, N2 = Nucleus.Table.ge(Nucleus.Table.rea(tsv2))

    bo_::Vector{Bool} = N1[1, :]

    N2 = N2[:, indexin(s1_, s3_)]

    update!(N2, standard_deviation)

    fu = if metric == "mean-difference"

        Nucleus.PairMetric.make_mean_difference

    elseif metric == "log-ratio"

        Nucleus.PairMetric.make_log_ratio

    elseif metric == "signal-to-noise-ratio"

        Nucleus.PairMetric.make_signal_to_noise_ratio

    end

    nu_ = map(nu_ -> Nucleus.PairMap.make(fu, bo_, nu_), eachrow(N2))

    Nucleus.Table.writ(
        joinpath(directory, "metric.tsv"),
        Nucleus.Table.make(st, s2_, [metric], reshape(nu_, :, 1)),
    )

    al = Algorithm.make(algorithm)

    di = Nucleus.Dictionary.rea(json)

    st__ = collect(values(di))

    ke_ = (mi = minimum, ma = maximum, pr = fraction)

    Result.writ(
        directory,
        al,
        s2_,
        nu_,
        collect(keys(di)),
        st__,
        Interface.make(al, s2_, nu_, st__; ke_...),
        if permutation == "set"

            Rando.make(number_of_permutations, seed, al, s2_, nu_, st__; ke_...)

        elseif permutation == "sample"

            Rando.make(number_of_permutations, seed, al, s2_, fu, bo_, N2, st__; ke_...)

        end,
        number_of_plots,
        split(more_plots, ';'),
        low,
        high,
    )

end

@main

end
