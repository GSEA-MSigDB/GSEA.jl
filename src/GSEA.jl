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

        n1 += n2 = Public.number_divergence(-, p2, p4, p3, p5)

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

        n1 += n3 = Public.number_divergence(-, p3, p5, p4, p6)

        if !isnothing(n2_)

            n2_[nd] = n3

        end

    end

    n1 / length(n1_)

end

########################################

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
                Public.number_divergence(-, r2, r1, r3, r3) -
                Public.number_divergence(-, l2, l1, l3, l3)

        if !isnothing(n2_)

            n2_[nd] = n3

        end

    end

    n1 / length(n1_)

end

########################################

function make_sort(a1_, n1_)

    in_ = findall(!isnan, n1_)

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
                    "line" => Dict("width" => 2, "color" => Public.TU),
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
                    # TODO
                    #"spikecolor" => Public.DA,
                ),
            ),
            d1,
        ),
    )

end

########################################
# TODO: Pick up

function write_enrichment(
    pa,
    al,
    s1_,
    s2_,
    N,
    s3_,
    st__,
    E,
    di = Dict{String, Any}();
    um = 2,
    ke_...,
)

    Public.write_heat(
        pa,
        s3_,
        s2_,
        E,
        Public.pair_merge(
            Dict(
                "title" => Dict("text" => "Enrichment"),
                "yaxis" => Dict("title" => Dict("text" => "Set")),
                "xaxis" => Dict("title" => Dict("text" => "Sample")),
            ),
            di,
        ),
    )

    in_ = findall(nu_ -> all(!isnan, nu_), eachrow(E))

    s3_ = s3_[in_]

    st__ = st__[in_]

    E = E[in_, :]

    for in_ in CartesianIndices(E)[Public.index_extreme(vec(E), um)]

        i3, i2 = Tuple(in_)

        st = s3_[i3]

        write_enrichment(
            "$(rsplit(pa, '.'; limit = 2)[1]).$(Public.text_2(E[in_])).$(s2_[i2]).$st.html",
            al,
            s1_,
            N[:, i2],
            st__[i3],
            Dict("title" => Dict("text" => st));
            ke_...,
        )

    end

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

function read_pair(js)

    di = Public.read_pair(js)

    collect(keys(di)), collect(values(di))

end

function number_z!(N, st)

    if iszero(st)

        return

    end

    for nu_ in eachcol(N)

        clamp!(Public.number_z!(nu_), -st, st)

    end

end

########################################

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `directory`:
  - `tsv`:
  - `json`:

# Options

  - `--standard-deviation`: For column-wise normalization. 0 skips normalization.
  - `--exponent`:
  - `--algorithm`: "S0" | "S0a" | "D2" | "D2f" | "D0f2f".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
  - `--number-of-plots`:
"""
@cast function data_rank(
    directory,
    tsv,
    json;
    standard_deviation::Float64 = 0.0,
    exponent::Float64 = 1.0,
    algorithm = "S0",
    minimum::Int = 1,
    maximum::Int = 1000,
    fraction::Float64 = 0.0,
    number_of_plots::Int = 2,
)

    al = make_algorithm(algorithm)

    st, s1_, s2_, N = Public.make_part(Public.read_table(tsv))

    number_z!(N, standard_deviation)

    s3_, st__ = read_pair(json)

    E = reduce(
        hcat,
        number_enrichment(
            al,
            s1_,
            nu_,
            st__;
            u1 = minimum,
            u2 = maximum,
            pr = fraction,
        ) for nu_ in eachcol(N)
    )

    fi = joinpath(directory, "result")

    Public.write_table("$fi.tsv", Public.make_table("Set", s3_, s2_, E))

    write_enrichment(
        "$fi.html",
        al,
        s1_,
        s2_,
        N,
        s3_,
        st__,
        E;
        um = number_of_plots,
    )

end

########################################

function number_random(um, se, al, s1_, nu_, st__; ke_...)

    R = Matrix{Float64}(undef, length(st__), um)

    um_ = map(s2_ -> length(intersect(s1_, s2_)), st__)

    seed!(se)

    @showprogress for nd in 1:um

        R[:, nd] = number_enrichment(
            al,
            s1_,
            nu_,
            map(um -> sample(s1_, um; replace = false), um_);
            ke_...,
        )

    end

    R

end

function number_random(um, se, al, st_, fu, bo_, N, st__; ke_...)

    R = Matrix{Float64}(undef, length(st__), um)

    seed!(se)

    @showprogress for nd in 1:um

        R[:, nd] = number_enrichment(
            al,
            st_,
            map(nu_ -> Public.make_2(fu, shuffle!(bo_), nu_), eachrow(N)),
            st__;
            ke_...,
        )

    end

    R

end

########################################

function number_normalization(n1, n2, n3)

    n1 / if n1 < 0

        -n2

    else

        n3

    end

end

function write_result(di, al, s1_, nu_, s2_, st__, n2_, R, um, s3_)

    N = Matrix{Float64}(undef, length(s2_), 4)

    N[:, 1] = n2_

    for i1 in axes(R, 1)

        m1, m2 = (mean(nu_) for nu_ in Public.number_sign(R[i1, :]))

        n2_[i1] = number_normalization(n2_[i1], m1, m2)

        for i2 in axes(R, 2)

            R[i1, i2] = number_normalization(R[i1, i2], m1, m2)

        end

    end

    N[:, 2] = n2_

    in_, pv_, qv_ = Public.number_significance(n2_, R)

    N[in_, 3] = pv_

    N[in_, 4] = qv_

    Public.write_table(
        joinpath(di, "result.tsv"),
        Public.make_table(
            "Set",
            s2_,
            ["Enrichment", "Normalized Enrichment", "P-Value", "Q-Value"],
            N,
        ),
    )

    # TODO: Use normalized enrichment
    in_ = findall(!isnan, N[:, 1])

    s2_ = s2_[in_]

    st__ = st__[in_]

    N = N[in_, :]

    for nd in unique!(
        vcat(
            Public.index_extreme(N[:, 1], um),
            filter!(!isnothing, indexin(s3_, s2_)),
        ),
    )

        st = s2_[nd]

        write_enrichment(
            joinpath(di, "$(Public.text_2(N[nd, 1])).$st.html"),
            al,
            s1_,
            nu_,
            st__[nd],
            Dict(Public.pair_title(st)),
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

  - `--exponent`:
  - `--algorithm`: "S0" | "S0a" | "D2" | "D2f" | "D0f2f".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
  - `--number-of-permutations`:
  - `--seed`:
  - `--number-of-plots`:
  - `--more-plots`: ;-separated set names.
"""
@cast function user_rank(
    directory,
    tsv,
    json;
    exponent::Float64 = 1.0,
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

    write_result(
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

  - `--standard-deviation`: For column-wise normalization. 0 skips normalization.
  - `--metric`: "signal-to-noise-ratio" | "mean-difference" | "log-ratio".
  - `--exponent`:
  - `--algorithm`: "S0" | "S0a" | "D2" | "D2f" | "D0f2f".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
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
    standard_deviation::Float64 = 0.0,
    metric = "signal-to-noise-ratio",
    exponent::Float64 = 1.0,
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

    _, _, s1_, N = Public.make_part(Public.read_table(tsv1))

    bo_::Vector{Bool} = view(N, 1, :)

    st, s2_, s3_, N = Public.make_part(Public.read_table(tsv2))

    N = N[:, indexin(s1_, s3_)]

    number_z!(N, standard_deviation)

    fu = if metric == "mean-difference"

        Public.number_difference

    elseif metric == "log-ratio"

        Public.number_ratio

    elseif metric == "signal-to-noise-ratio"

        Public.number_signal

    end

    # TODO: Consider Match

    nu_ = map(nu_ -> Public.make_2(fu, bo_, nu_), eachrow(N))

    Public.write_table(
        joinpath(directory, "metric.tsv"),
        Public.make_table(st, s2_, [metric], reshape(nu_, :, 1)),
    )

    al = make_algorithm(algorithm)

    s3_, st__ = read_pair(json)

    ke_ = (u1 = minimum, u2 = maximum, pr = fraction)

    write_result(
        directory,
        al,
        s2_,
        nu_,
        s3_,
        st__,
        number_enrichment(al, s2_, nu_, st__; ke_...),
        if permutation == "set"

            number_random(
                number_of_permutations,
                seed,
                al,
                s2_,
                nu_,
                st__;
                ke_...,
            )

        elseif permutation == "sample"

            number_random(
                number_of_permutations,
                seed,
                al,
                s2_,
                fu,
                bo_,
                N,
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
