module GSEA

const P1 = pkgdir(GSEA, "in")

const P2 = pkgdir(GSEA, "ou")

# ------------------------------------ #

using Comonicon: @cast, @main

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, sample

using Public

########################################

struct S0 end

struct S0a end

struct D2 end

struct D2f end

struct D0f2f end

########################################

function number_delta(::Union{S0, S0a}, nu_, bo_)

    n1 = n2 = 0.0

    for nd in eachindex(nu_)

        if bo_[nd]

            n2 += abs(nu_[nd])

        else

            n1 += 1.0

        end

    end

    -inv(n1), inv(n2)

end

function number_delta(::Union{D2, D2f, D0f2f}, nu_, bo_)

    n1 = n2 = 0.0

    for nd in eachindex(nu_)

        n2 += n3 = abs(nu_[nd])

        if bo_[nd]

            n1 += n3

        end

    end

    inv(n1), inv(n2)

end

function number_delta(p1, p2)

    inv(inv(p2) - inv(p1))

end

########################################

function number_enrichment!(al::S0, n1_, bo_, n2_ = nothing)

    p1, p2 = number_delta(al, n1_, bo_)

    n1 = n2 = n3 = 0.0

    for nd in eachindex(n1_)

        n1 += if bo_[nd]

            p2 * abs(n1_[nd])

        else

            p1

        end

        if !isnothing(n2_)

            n2_[nd] = n1

        end

        n4 = abs(n1)

        if n2 < n4

            n2 = n4

            n3 = n1

        end

    end

    n3

end

function number_enrichment!(al::S0a, n1_, bo_, n2_ = nothing)

    p1, p2 = number_delta(al, n1_, bo_)

    n1 = n2 = 0.0

    for nd in eachindex(n1_)

        n2 += n1 += if bo_[nd]

            p2 * abs(n1_[nd])

        else

            p1

        end

        if !isnothing(n2_)

            n2_[nd] = n1

        end

    end

    n2 / length(n1_)

end

########################################

function number_eps(nu)

    max(eps(), nu)

end

########################################

function number_enrichment!(al::D2, n1_, bo_, n2_ = nothing)

    um = length(n1_)

    p1, _ = number_delta(al, n1_, bo_)

    p2 = p3 = eps()

    p4 = p5 = 1.0

    p6 = n1 = 0.0

    p5 += p7 = inv(um)

    for nd in eachindex(n1_)

        p4 = number_eps(p4 - p6)

        p5 = number_eps(p5 - p7)

        p2 += p6 = if bo_[nd]

            p1 * abs(n1_[nd])

        else

            0.0

        end

        p3 += p7

        n1 += n2 = Public.number_divergence(-, p2, p3, p4, p5)

        if !isnothing(n2_)

            n2_[nd] = n2

        end

    end

    n1 / um

end

function number_enrichment!(al::D2f, n1_, bo_, n2_ = nothing)

    p1, p2 = number_delta(al, n1_, bo_)

    p3 = p4 = eps()

    p5 = p6 = 1.0

    p7 = p8 = n1 = 0.0

    for nd in eachindex(n1_)

        p5 = number_eps(p5 - p7)

        p6 = number_eps(p6 - p8)

        n2 = abs(n1_[nd])

        p3 += p7 = if bo_[nd]

            p1 * n2

        else

            0.0

        end

        p4 += p8 = p2 * n2

        n1 += n3 = Public.number_divergence(-, p3, p4, p5, p6)

        if !isnothing(n2_)

            n2_[nd] = n3

        end

    end

    n1 / length(n1_)

end

function number_enrichment!(al::D0f2f, n1_, bo_, n2_ = nothing)

    p1, p2 = number_delta(al, n1_, bo_)

    p3 = number_delta(p1, p2)

    r1 = r2 = r3 = eps()

    l1 = l2 = l3 = 1.0

    c1 = c2 = c3 = n1 = 0.0

    for nd in eachindex(n1_)

        l1 = number_eps(l1 - c1)

        l2 = number_eps(l2 - c2)

        l3 = number_eps(l3 - c3)

        n2 = abs(n1_[nd])

        c1, c2 = if bo_[nd]

            0.0, p1 * n2

        else

            p3 * n2, 0.0

        end

        r1 += c1

        r2 += c2

        r3 += c3 = p2 * n2

        n1 +=
            n3 =
                Public.number_divergence(-, r2, r3, r1, r3) -
                Public.number_divergence(-, l2, l3, l1, l3)

        if !isnothing(n2_)

            n2_[nd] = n3

        end

    end

    n1 / length(n1_)

end

########################################

function make_sort(a1_, n1_)

    in_ = findall(isfinite, n1_)

    a2_ = a1_[in_]

    n2_ = n1_[in_]

    # TODO: Make increasing
    sortperm!(in_, n2_; rev = true)

    a2_[in_], n2_[in_]

end

########################################

function number_enrichment(al, s1_, n1_, st__; u1 = 1, u2 = 1000, pr = 0)

    s2_, n2_ = make_sort(s1_, n1_)

    u3 = length(s2_)

    bo_ = falses(u3)

    di = Dict(s2_[nd] => nd for nd in 1:u3)

    u4 = length(st__)

    n3_ = Vector{Float64}(undef, u4)

    for nd in 1:u4

        s3_ = st__[nd]

        for st in s3_

            if haskey(di, st)

                bo_[di[st]] = true

            end

        end

        u5 = count(bo_)

        n3_[nd] = if u1 <= u5 <= u2 && pr <= u5 / length(s3_)

            number_enrichment!(al, n2_, bo_)

        else

            NaN

        end

        fill!(bo_, false)

    end

    n3_

end

########################################

function write_enrichment(pa, al, s1_, n1_, s2_, d1 = Dict{String, Any}())

    s3_, n2_ = make_sort(s1_, n1_)

    um = length(s3_)

    n3_ = Vector{Float64}(undef, um)

    bo_ = map(in(s2_), s3_)

    st = Public.text_4(number_enrichment!(al, n2_, bo_, n3_))

    d2 = Dict(
        "mode" => "lines",
        "line" => Dict("width" => 0),
        "fill" => "tozeroy",
    )

    i1_ = findall(<(0), n2_)

    i2_ = findall(>=(0), n2_)

    Public.write_plotly(
        pa,
        (
            merge(
                d2,
                Dict(
                    "yaxis" => "y3",
                    "y" => n3_,
                    "x" => s3_,
                    "fillcolor" => "#07fa07",
                ),
            ),
            Dict(
                "yaxis" => "y2",
                "y" => zeros(count(bo_)),
                "x" => s3_[bo_],
                "mode" => "markers",
                "marker" => Dict(
                    "symbol" => "line-ns",
                    "size" => 24,
                    "line" => Dict(
                        "width" => 2,
                        "color" => Public.text_color(Public.IN, 0.8),
                    ),
                ),
            ),
            merge(
                d2,
                Dict(
                    "y" => n2_[i1_],
                    "x" => s3_[i1_],
                    "fillcolor" => Public.BL,
                ),
            ),
            merge(
                d2,
                Dict(
                    "y" => n2_[i2_],
                    "x" => s3_[i2_],
                    "fillcolor" => Public.RE,
                ),
            ),
        ),
        Public.pair_merge(
            Dict(
                "showlegend" => false,
                Public.pair_title("", "Enrichment = <b>$st</b>"),
                "yaxis3" => Dict(
                    "domain" => (0.328, 1),
                    Public.pair_title("Δ Enrichment"),
                ),
                "yaxis2" => Dict(
                    "domain" => (0.248, 0.32),
                    Public.pair_title("Set"),
                    "tickvals" => (),
                ),
                "yaxis" =>
                    Dict("domain" => (0, 0.24), Public.pair_title("Score")),
                "xaxis" => Dict(
                    Public.pair_title("Feature ($um)"),
                    "showspikes" => true,
                    "spikemode" => "across",
                    "spikedash" => "solid",
                    "spikethickness" => -1,
                    "spikecolor" => Public.DA,
                ),
            ),
            d1,
        ),
    )

end

########################################

function make_algorithm(st)

    if st == "S0"

        S0()

    elseif st == "S0a"

        S0a()

    elseif st == "D2"

        D2()

    elseif st == "D2f"

        D2f()

    elseif st == "D0f2f"

        D0f2f()

    end

end

function read_pair(pa)

    di::Dict{String, Vector{String}} = Public.read_pair(pa)

    st_ = collect(keys(di))

    in_ = sortperm(st_)

    st_[in_], collect(values(di))[in_]

end

########################################

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `directory`:
  - `tsv`:
  - `json`:

# Options

  - `--algorithm`: "S0" | "S0a" | "D2" | "D2f" | "D0f2f".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members.
  - `--number-of-plots`:
"""
@cast function data_rank(
    directory,
    tsv,
    json;
    algorithm = "S0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Float64 = 0.0,
    number_of_plots::Int = 2,
)

    al = make_algorithm(algorithm)

    _, s1_, s2_, N1 = Public.make_part(Public.read_table(tsv))

    s3_, s1__ = read_pair(json)

    N2 = reduce(
        hcat,
        number_enrichment(
            al,
            s1_,
            nu_,
            s1__;
            u1 = minimum,
            u2 = maximum,
            pr = fraction,
        ) for nu_ in eachcol(N1)
    )

    pa = joinpath(directory, "enrichment")

    Public.write_table("$pa.tsv", Public.make_table("Set", s3_, s2_, N2))

    Public.write_heat(
        "$pa.html",
        s3_,
        s2_,
        N2,
        Dict(
            Public.pair_title("Enrichment"),
            "yaxis" => Dict(Public.pair_title("Set")),
            "xaxis" => Dict(Public.pair_title("Sample")),
        ),
    )

    bo_ = map(nu_ -> any(isfinite, nu_), eachrow(N2))

    s4_ = s3_[bo_]

    s2__ = s1__[bo_]

    N3 = N2[bo_, :]

    for in_ in
        CartesianIndices(N3)[Public.index_extreme(vec(N3), number_of_plots)]

        i1, i2 = Tuple(in_)

        s1 = s4_[i1]

        s2 = s2_[i2]

        s3 = Public.text_2(N3[in_])

        write_enrichment(
            "$pa.$s1.$s2.$s3.html",
            al,
            s1_,
            N1[:, i2],
            s2__[i1],
            Dict(Public.pair_title(s1)),
        )

    end

end

########################################

function number_random(um, nu, al, s1_, nu_, st__; ke_...)

    N = Matrix{Float64}(undef, length(st__), um)

    um_ = map(s2_ -> length(intersect(s1_, s2_)), st__)

    seed!(nu)

    @showprogress for nd in 1:um

        N[:, nd] = number_enrichment(
            al,
            s1_,
            nu_,
            map(um -> sample(s1_, um; replace = false), um_);
            ke_...,
        )

    end

    N

end

function number_random!(um, nu, al, st_, fu, bo_, N1, st__; ke_...)

    N2 = Matrix{Float64}(undef, length(st__), um)

    seed!(nu)

    @showprogress for nd in 1:um

        N2[:, nd] = number_enrichment(
            al,
            st_,
            map(nu_ -> Public.make_2(fu, shuffle!(bo_), nu_), eachrow(N1)),
            st__;
            ke_...,
        )

    end

    N2

end

########################################

function number_normalization(n1, n2, n3)

    n1 / if n1 < 0

        -n2

    else

        n3

    end

end

function write_table(di, al, s1_, n1_, s2_, s1__, n2_, N1, u1, s3_)

    u2 = length(s2_)

    N2 = Matrix{Float64}(undef, u2, 4)

    N2[:, 1] = n2_

    for i1 in 1:u2

        n1, n2 = (mean(n6_) for n6_ in Public.number_sign(N1[i1, :]))

        n2_[i1] = number_normalization(n2_[i1], n1, n2)

        for i2 in axes(N1, 2)

            N1[i1, i2] = number_normalization(N1[i1, i2], n1, n2)

        end

    end

    N2[:, 2] = n2_

    in_, n3_, n4_ = Public.number_significance(n2_, N1)

    N2[in_, 3] = n3_

    N2[in_, 4] = n4_

    Public.write_table(
        joinpath(di, "enrichment.tsv"),
        Public.make_table(
            "Set",
            s2_,
            ["Enrichment", "Normalized enrichment", "P value", "Q value"],
            N2,
        ),
    )

    bo_ = map(isfinite, n2_)

    s4_ = s2_[bo_]

    s2__ = s1__[bo_]

    n5_ = n2_[bo_]

    for nd in unique!(
        vcat(
            Public.index_extreme(n5_, u1),
            filter!(!isnothing, indexin(s3_, s4_)),
        ),
    )

        s1 = s4_[nd]

        s2 = Public.text_2(n5_[nd])

        write_enrichment(
            joinpath(di, "$s2.$s1.html"),
            al,
            s1_,
            n1_,
            s2__[nd],
            Dict(Public.pair_title(s1)),
        )

    end

end

########################################

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `directory`:
  - `tsv`:
  - `json`:

# Options

  - `--algorithm`: "S0" | "S0a" | "D2" | "D2f" | "D0f2f".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members.
  - `--number-of-permutations`:
  - `--seed`:
  - `--number-of-plots`:
  - `--more-plots`: ;-separated set names.
"""
@cast function user_rank(
    directory,
    tsv,
    json;
    algorithm = "S0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Float64 = 0.0,
    number_of_permutations::Int = 100,
    seed::Int = 20150603,
    number_of_plots::Int = 2,
    more_plots = "",
)

    s1_, nu_ = eachcol(Public.read_table(tsv)[!, 1:2])

    al = make_algorithm(algorithm)

    s2_, st__ = read_pair(json)

    ke_ = (u1 = minimum, u2 = maximum, pr = fraction)

    write_table(
        directory,
        al,
        s1_,
        nu_,
        s2_,
        st__,
        number_enrichment(al, s1_, nu_, st__; ke_...),
        number_random(number_of_permutations, seed, al, s1_, nu_, st__; ke_...),
        number_of_plots,
        split(more_plots, ';'),
    )

end

########################################

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `directory`:
  - `tsv1`:
  - `tsv2`:
  - `json`:

# Options

  - `--metric`: "signal-to-noise-ratio" | "mean-difference" | "log-ratio".
  - `--algorithm`: "S0" | "S0a" | "D2" | "D2f" | "D0f2f".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members.
  - `--permutation`: "sample" | "set".
  - `--number-of-permutations`:
  - `--seed`:
  - `--number-of-plots`:
  - `--more-plots`: ;-separated set names.
"""
@cast function metric_rank(
    directory,
    tsv1,
    tsv2,
    json;
    metric = "signal-to-noise-ratio",
    algorithm = "S0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Float64 = 0.0,
    permutation = "sample",
    number_of_permutations::Int = 100,
    seed::Int = 20150603,
    number_of_plots::Int = 2,
    more_plots = "",
)

    _, _, s1_, N1 = Public.make_part(Public.read_table(tsv1))

    bo_::BitVector = N1[1, :]

    st, s2_, s3_, N1 = Public.make_part(Public.read_table(tsv2))

    N2 = N1[:, indexin(s1_, s3_)]

    fu = if metric == "mean-difference"

        Public.number_difference

    elseif metric == "log-ratio"

        Public.number_ratio

    elseif metric == "signal-to-noise-ratio"

        Public.number_signal

    end

    n1_ = map(n2_ -> Public.make_2(fu, bo_, n2_), eachrow(N2))

    Public.write_table(
        joinpath(directory, "metric.tsv"),
        Public.make_table(st, s2_, [metric], reshape(n1_, :, 1)),
    )

    al = make_algorithm(algorithm)

    s4_, st__ = read_pair(json)

    ke_ = (u1 = minimum, u2 = maximum, pr = fraction)

    write_table(
        directory,
        al,
        s2_,
        n1_,
        s4_,
        st__,
        number_enrichment(al, s2_, n1_, st__; ke_...),
        if permutation == "set"

            number_random(
                number_of_permutations,
                seed,
                al,
                s2_,
                n1_,
                st__;
                ke_...,
            )

        elseif permutation == "sample"

            number_random!(
                number_of_permutations,
                seed,
                al,
                s2_,
                fu,
                bo_,
                N2,
                st__;
                ke_...,
            )

        end,
        number_of_plots,
        split(more_plots, ';'),
    )

end

@main

end
