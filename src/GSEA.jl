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
# TODO: Rename

struct KS0 end

struct A0 end

struct DA2 end

struct DA2W end

struct DA2W0W end

########################################

function number_delta(::Union{KS0, A0}, nu_, bo_)

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

function number_delta(::Union{DA2, DA2W, DA2W0W}, nu_, bo_)

    n1 = n2 = 0.0

    for nd in eachindex(nu_)

        n2 += ab = abs(nu_[nd])

        if bo_[nd]

            n1 += ab

        end

    end

    inv(n1), inv(n2)

end

function number_delta(p1, p2)

    inv(inv(p2) - inv(p1))

end

########################################

function number_enrichment!(al::KS0, n1_, bo_, n2_ = nothing)

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

function number_enrichment!(al::A0, n1_, bo_, n2_ = nothing)

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

function number_enrichment!(al::DA2, n1_, bo_, n2_ = nothing)

    pr ,_ = number_delta(al, n1_, bo_)

    n1 = n2 = eps()

    n3 = n4 = 1.0

    n5 = n7 = 0.0

    n4 += n6 = inv(length(n1_))

    for nd in eachindex(n1_)

        n3 = number_eps(n3 - n5)

        n4 = number_eps(n4 - n6)

        n1 += n5 = if bo_[nd] 

            pr * abs(n1_[nd]) 

            else

                0.0

            end

        n2 += n6

        n7 += n8 = Public.number_divergence(-, n1, n3, n2, n4)

        if !isnothing(n2_)

            n2_[nd] = n8

        end

    end

    n7 / length(n1_)

end

function number_enrichment!(al::DA2W, n1_, bo_, n2_ = nothing)

    p1, p2 = number_delta(al, n1_, bo_)

    n1 = n2 = eps()

    n3 = n4 = 1.0

    n5 = n6 = n7 = 0.0

    for nd in eachindex(n1_)

        n8 = abs(n1_[nd])

        n3 = number_eps(n3 - n5)

        n4 = number_eps(n4 - n6)

        n1 += n5 = if bo_[nd] 

            p1 * n8 

            else

                0.0

            end

        n2 += n6 = p2 * n8

        n7 += n9 = Public.number_divergence(-, n1, n3, n2, n4)

        if !isnothing(n2_)

            n2_[nd] = n9

        end

    end

    n7 / length(n1_)

end

########################################

function number_enrichment!(al::DA2W0W, n1_, bo_, n2_ = nothing)

    p1, p2 = number_delta(al, n1_, bo_)

    p3 = number_delta(p1, p2)

    r0 = r1 = r2 = eps()

    l0 = l1 = l2 = 1.0

    e0 = e1 = e2 = su = 0.0

    for nd in eachindex(n1_)

        ab = abs(n1_[nd])

        l0 = number_eps(l0 - e0)

        l1 = number_eps(l1 - e1)

        l2 = number_eps(l2 - e2)

        r0 += e0 = bo_[nd] ? 0.0 : p3 * ab

        r1 += e1 = bo_[nd] ? p1 * ab : 0.0

        r2 += e2 = p2 * ab

        su +=
            cu =
                Public.number_divergence(-, r1, r0, r2, r2) -
                Public.number_divergence(-, l1, l0, l2, l2)

        if !isnothing(n2_)

            n2_[nd] = cu

        end

    end

    su / length(n1_)

end

########################################

function make_sort(a1_, n1_)

    in_ = findall(!isnan, n1_)

    a2_ = a1_[in_]

    n2_ = n1_[in_]

    sortperm!(in_, n2_; rev = true)

    a2_[in_], n2_[in_]

end

########################################

function number_enrichment(al, s1_, nu_, st__; mi = 1, ma = 1000, pr = 0)

    s1_, nu_ = make_sort(s1_, nu_)

    di = Dict(s1_[nd] => nd for nd in eachindex(s1_))

    bo_ = falses(length(s1_))

    en_ = Vector{Float64}(undef, length(st__))

    for i1 in eachindex(st__)

        s2_ = st__[i1]

        for x1 in s2_

            i2 = get(di, x1, nothing)

            if !isnothing(i2)

                bo_[i2] = true

            end

        end

        um = sum(bo_)

        en_[i1] =
            um < mi || ma < um || um / length(s2_) < pr ? NaN :
            number_enrichment!(al, nu_, bo_)

        fill!(bo_, false)

    end

    en_

end

########################################

function write_enrichment(
    ht,
    al,
    s1_,
    nu_,
    s2_,
    la = Dict{String, Any}();
    t1 = "Score",
    t2 = "Low",
    t3 = "High",
)

    s1_, nu_ = make_sort(s1_, nu_)

    um = length(s1_)

    tr = Dict(
        "mode" => "lines",
        "line" => Dict("width" => 0),
        "fill" => "tozeroy",
    )

    xc_ = 1:um

    i1_ = findall(<(0), nu_)

    i2_ = findall(>=(0), nu_)

    bo_ = map(in(s2_), s1_)

    cu_ = Vector{Float64}(undef, um)

    en = number_enrichment!(al, nu_, bo_, cu_)

    an = Dict(
        "y" => 0,
        "font" => Dict("size" => 16),
        "borderpad" => 4.8,
        "borderwidth" => 2.64,
        "bordercolor" => Public.LI,
        "showarrow" => false,
    )

    po = um * 0.008

    Public.write_plotly(
        ht,
        (
            merge(
                tr,
                Dict(
                    "y" => nu_[i1_],
                    "x" => xc_[i1_],
                    "text" => s1_[i1_],
                    "fillcolor" => Public.BL,
                ),
            ),
            merge(
                tr,
                Dict(
                    "y" => nu_[i2_],
                    "x" => xc_[i2_],
                    "text" => s1_[i2_],
                    "fillcolor" => Public.RE,
                ),
            ),
            Dict(
                "yaxis" => "y2",
                "y" => zeros(sum(bo_)),
                "x" => xc_[bo_],
                "mode" => "markers",
                "marker" => Dict(
                    "symbol" => "line-ns",
                    "size" => 24,
                    "line" => Dict("width" => 2, "color" => "#000000cc"),
                ),
                "text" => s1_[bo_],
                "hoverinfo" => "x+text",
            ),
            merge(
                tr,
                Dict(
                    "yaxis" => "y3",
                    "y" => cu_,
                    "x" => xc_,
                    "text" => s1_,
                    "fillcolor" => "#07fa07",
                ),
            ),
        ),
        Public.pair_merge(
            Dict(
                "showlegend" => false,
                "yaxis" =>
                    Dict("domain" => (0, 0.24), "title" => Dict("text" => t1)),
                "yaxis2" => Dict(
                    "domain" => (0.248, 0.32),
                    "title" => Dict("text" => "Set"),
                    "tickvals" => (),
                ),
                "yaxis3" => Dict(
                    "domain" => (0.328, 1),
                    "title" => Dict("text" => "Δ Enrichment"),
                ),
                "xaxis" => Dict(
                    "title" => Dict("text" => "Feature ($um)"),
                    "showspikes" => true,
                    "spikemode" => "across",
                    "spikedash" => "solid",
                    "spikethickness" => -1,
                    "spikecolor" => "#000000",
                ),
                "annotations" => (
                    Dict(
                        "yref" => "paper",
                        "xref" => "paper",
                        "y" => 1.056,
                        "text" => "Enrichment = <b>$(Public.text_4(en))</b>",
                        "font" => Dict("size" => 24, "color" => "#000000"),
                        "showarrow" => false,
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => 1 - po,
                            "xanchor" => "right",
                            "text" => t3,
                            "font" => Dict("color" => Public.RE),
                        ),
                    ),
                    merge(
                        an,
                        Dict(
                            "x" => um + po,
                            "xanchor" => "left",
                            "text" => t2,
                            "font" => Dict("color" => Public.BL),
                        ),
                    ),
                ),
            ),
            la,
        ),
    )

end

########################################

function write_enrichment(
    ht,
    al,
    s1_,
    s2_,
    N,
    s3_,
    st__,
    E,
    la = Dict{String, Any}();
    um = 2,
    ke_...,
)

    Public.write_heat(
        ht,
        s3_,
        s2_,
        E,
        Public.pair_merge(
            Dict(
                "title" => Dict("text" => "Enrichment"),
                "yaxis" => Dict("title" => Dict("text" => "Set")),
                "xaxis" => Dict("title" => Dict("text" => "Sample")),
            ),
            la,
        ),
    )

    in_ = findall(en_ -> all(!isnan, en_), eachrow(E))

    s3_ = s3_[in_]

    st__ = st__[in_]

    E = E[in_, :]

    for in_ in CartesianIndices(E)[Public.index_extreme(vec(E), um)]

        i3, i2 = Tuple(in_)

        st = s3_[i3]

        write_enrichment(
            "$(rsplit(ht, '.'; limit = 2)[1]).$(Public.text_2(E[in_])).$(s2_[i2]).$st.html",
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

# #numeric
# #Score
# 1 2 4 8 1
#
# 6 2 1
# #Aa Bb
# Aa Aa Aa Bb Bb Bb
function read_cls(pa)

    s1, s2, s3 = readlines(pa)

    s4 = s2[2:end]

    s1_ = split(s3)

    s5, nu_ = if s1 == "#numeric"

        s4, map(s6 -> parse(Float64, s6), s1_)

    else

        s2_ = split(s4)

        di = Dict(s2_[nd] => nd for nd in eachindex(s2_))

        join(s2_, " vs "), map(s6 -> di[s6], s1_)

    end

    "Phenotype",
    s5,
    map!(nd -> "Sample$nd", s1_, eachindex(s1_)),
    reshape(nu_, 1, :)

end

function read_gct(pa)

    Public.make_part(Public.read_table(pa; header = 3, drop = ["Description"]))

end

function read_gmt(pa)

    di = Dict{String, Vector{String}}()

    for s1 in eachline(pa)

        st_ = split(s1, '\t')

        di[st_[1]] = filter!(!isempty, st_[3:end])

    end

    di

end

########################################

function make_algorithm(st)

    if st == "ks0"

        KS0()

    elseif st == "a0"

        A0()

    elseif st == "da2"

        DA2()

    elseif st == "da2w"

        DA2W()

    elseif st == "da2w0w"

        DA2W0W()

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
  - `--algorithm`: "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
  - `--number-of-plots`:
  - `--score`:
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
    score = "Data",
    low = "Low",
    high = "High",
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
            mi = minimum,
            ma = maximum,
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
        t1 = score,
        t2 = low,
        t3 = high,
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

function write_result(di, al, s1_, nu_, s2_, st__, en_, R, um, s3_, t1, t2, t3)

    N = Matrix{Float64}(undef, length(s2_), 4)

    N[:, 1] = en_

    for i1 in axes(R, 1)

        m1, m2 = (mean(nu_) for nu_ in Public.number_sign(R[i1, :]))

        en_[i1] = number_normalization(en_[i1], m1, m2)

        for i2 in axes(R, 2)

            R[i1, i2] = number_normalization(R[i1, i2], m1, m2)

        end

    end

    N[:, 2] = en_

    in_, pv_, qv_ = Public.number_significance(en_, R)

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
            Dict("title" => Dict("text" => st));
            t1,
            t2,
            t3,
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
  - `--algorithm`: "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
  - `--number-of-permutations`:
  - `--seed`:
  - `--number-of-plots`:
  - `--more-plots`: ;-separated set names.
  - `--score`:
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
    score = "",
    low = "Low",
    high = "High",
)

    s1_, nu_ = eachcol(Public.read_table(tsv)[!, 1:2])

    al = make_algorithm(algorithm)

    s2_, st__ = read_pair(json)

    ke_ = (mi = minimum, ma = maximum, pr = fraction)

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
        score,
        low,
        high,
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
  - `--algorithm`: "ks0" | "a0" | "da2" | "da2w" | "da2w0w".
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum fraction of set members present.
  - `--permutation`: "sample" | "set".
  - `--number-of-permutations`:
  - `--seed`:
  - `--number-of-plots`:
  - `--more-plots`: ;-separated set names.
  - `--score`:
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
    score = metric,
    low = "Low",
    high = "High",
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

    ke_ = (mi = minimum, ma = maximum, pr = fraction)

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
        score,
        low,
        high,
    )

end

@main

end
