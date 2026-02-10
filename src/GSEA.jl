module GSEA

const P1 = pkgdir(GSEA, "in")

const P2 = pkgdir(GSEA, "ou")

# ------------------------------------ #

using CSV: read, write as write2

using Comonicon: @cast, @main

using DataFrames: DataFrame, insertcols!

using JSON: json, parsefile

using MultipleTesting: BenjaminiHochberg, adjust

using Printf: @sprintf

using ProgressMeter: @showprogress

using Random: seed!, shuffle!

using StatsBase: mean, mean_and_std, sample

########################################

function text_2(nu)

    @sprintf "%.2g" nu

end

function text_4(nu)

    @sprintf "%.4g" nu

end

########################################

function number_divergence(p1, p2)

    p1 * log2(p1 / p2)

end

function number_divergence(p1, p2, p3, p4, fu)

    fu(number_divergence(p1, p2), number_divergence(p3, p4))

end

########################################

function index_extreme(an_, u1)

    u2 = length(an_)

    sortperm(an_)[if u2 <= 2 * u1

        1:u2

    else

        vcat(1:u1, (u2 - u1 + 1):u2)

    end]

end

########################################

function number_significance(n1_, n2_, fu)

    u1 = length(n1_)

    u2 = length(n2_)

    if iszero(u1) || iszero(u2)

        return fill(NaN, u1), fill(NaN, u1)

    end

    pr_ = map(nu -> max(1, count(fu(nu), n2_)) / u2, n1_)

    pr_, adjust(pr_, BenjaminiHochberg())

end

# TODO: Use both signs
function number_significance(n1_, n2_)

    i1_ = findall(<(0), n1_)

    i2_ = findall(>=(0), n1_)

    n3_, n4_ = number_significance(n1_[i1_], filter(<(0), n2_), <=)

    n5_, n6_ = number_significance(n1_[i2_], filter(!<(0), n2_), >=)

    vcat(i1_, i2_), vcat(n3_, n5_), vcat(n4_, n6_)

end

function number_difference(n1_, n2_)

    mean(n2_) - mean(n1_)

end

function number_ratio(n1_, n2_)

    log2(mean(n2_) / mean(n1_))

end

function number_signal(n1_, n2_)

    n1, n2 = mean_and_std(n1_)

    n3, n4 = mean_and_std(n2_)

    (n3 - n1) / (max(0.2 * abs(n1), n2) + max(0.2 * abs(n3), n4))

end

########################################

function make_function(in_, an_, fu)

    fu(an_[findall(isone, in_)], an_[findall(==(2), in_)])

end

########################################

function read_open(pa)

    try

        run(`open --background $pa`)

    catch

        @info "⚠️ Failed to open $pa"

    end

end

########################################

function read_table(pa; ke_...)

    read(pa, DataFrame; ke_...)

end

function write_table(pa, D)

    write2(pa, D; delim = '\t')

end

function table_part(D)

    st_ = names(D)

    st_[1], D[!, 1], st_[2:end], Matrix(D[!, 2:end])

end

function table_part(st, s1_, s2_, A)

    insertcols!(DataFrame(A, s2_), 1, st => s1_)

end

########################################

function write_html(p1, pa_, s1, he = "#27221f")

    p2 = if isempty(p1)

        tempname(; cleanup = false, suffix = ".html")

    else

        p1

    end

    s2 = join("<script src=\"$p3\"></script>\n" for p3 in pa_)

    write(
        p2,
        """
        <!DOCTYPE html>
        <html lang="en">
        <head>
          <meta charset="utf-8">
        </head>
        $s2
        <body style="margin:0; background:$he; min-height:100vh; display:flex; justify-content:center; align-items:center">
          <div id="write_html" style="height:88vh; width:88vw"></div>
        </body>
        <script>
        $s1
        </script>
        </html>""",
    )

    read_open(p2)

end

########################################

function pair_font(nu)

    "font" => Dict("size" => nu)

end

function pair_title(st)

    "title" => Dict("text" => st)

end

function pair_title(s1, s2)

    "title" => Dict("text" => s1, "subtitle" => Dict("text" => s2))

end

function write_plotly(pa, di_, d1 = Dict(), d2 = Dict())

    s1 = json(di_; allownan = true)

    d3 = Dict(
        "automargin" => true,
        "title" => Dict(pair_font(24)),
        "zeroline" => false,
        "showgrid" => false,
    )

    s2 = json(
        merge(
            Dict(
                "template" => Dict(
                    "data" => Dict(
                        "scatter" => (Dict("cliponaxis" => false),),
                        "heatmap" => (
                            Dict(
                                "colorbar" => Dict(
                                    "lenmode" => "pixels",
                                    "len" => 240,
                                    "thickness" => 16,
                                    "outlinewidth" => 0,
                                ),
                            ),
                        ),
                    ),
                    "layout" => Dict(
                        "title" => Dict(pair_font(32)),
                        "yaxis" => d3,
                        "xaxis" => d3,
                        "legend" => Dict(pair_font(16)),
                    ),
                ),
            ),
            d1,
        ),
    )

    s3 = json(d2)

    write_html(
        pa,
        ("https://cdn.plot.ly/plotly-3.3.0.min.js",),
        "Plotly.newPlot(\"write_html\", $s1, $s2, $s3)",
    )

end

########################################

struct S0 end

struct S0a end

struct D2 end

struct D2w end

struct DD end

########################################

function number_weight(::Union{S0, S0a}, nu_, bo_)

    um = 0

    nu = 0.0

    for nd in eachindex(nu_)

        if bo_[nd]

            nu += abs(nu_[nd])

        else

            um += 1

        end

    end

    -inv(um), inv(nu)

end

function number_weight(::Union{D2, D2w, DD}, nu_, bo_)

    n1 = n2 = 0.0

    for nd in eachindex(nu_)

        n2 += n3 = abs(nu_[nd])

        if bo_[nd]

            n1 += n3

        end

    end

    inv(n1), inv(n2)

end

########################################

function number_enrichment!(al::S0, n1_, bo_, n2_ = nothing)

    p1, p2 = number_weight(al, n1_, bo_)

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

    p1, p2 = number_weight(al, n1_, bo_)

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

# TODO: Delete
function number_eps(nu)

    max(eps(), nu)

end

function number_enrichment!(al::D2, n1_, bo_, n2_ = nothing)

    um = length(n1_)

    p1, = number_weight(al, n1_, bo_)

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

        n1 += n2 = number_divergence(p2, p3, p4, p5, -)

        if !isnothing(n2_)

            n2_[nd] = n2

        end

    end

    n1 / um

end

function number_enrichment!(al::D2w, n1_, bo_, n2_ = nothing)

    p1, p2 = number_weight(al, n1_, bo_)

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

        n1 += n3 = number_divergence(p3, p4, p5, p6, -)

        if !isnothing(n2_)

            n2_[nd] = n3

        end

    end

    n1 / length(n1_)

end

function number_enrichment!(al::DD, n1_, bo_, n2_ = nothing)

    p1, p2 = number_weight(al, n1_, bo_)

    p3 = inv(inv(p2) - inv(p1))

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
                number_divergence(r2, r3, r1, r3, -) -
                number_divergence(l2, l3, l1, l3, -)

        if !isnothing(n2_)

            n2_[nd] = n3

        end

    end

    n1 / length(n1_)

end

########################################

function make_score(st_, n1_)

    bo_ = map(isfinite, n1_)

    n2_ = n1_[bo_]

    # TODO: Delete rev
    in_ = sortperm(n2_; rev = true)

    st_[bo_][in_], n2_[in_]

end

function write_enrichment(pa, al, s1_, n1_, s2_, s1 = "Mountain plot")

    s3_, n2_ = make_score(s1_, n1_)

    um = length(s3_)

    bo_ = map(in(s2_), s3_)

    n3_ = Vector{Float64}(undef, um)

    s2 = text_4(number_enrichment!(al, n2_, bo_, n3_))

    di = Dict(
        "mode" => "lines",
        "line" => Dict("width" => 0),
        "fill" => "tozeroy",
    )

    i1_ = findall(<(0), n2_)

    i2_ = findall(>=(0), n2_)

    write_plotly(
        pa,
        (
            merge(
                di,
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
                    "line" => Dict("width" => 2, "color" => "#4e40d8b8"),
                ),
            ),
            merge(
                di,
                Dict(
                    "y" => n2_[i1_],
                    "x" => s3_[i1_],
                    "fillcolor" => "#1992ff",
                ),
            ),
            merge(
                di,
                Dict(
                    "y" => n2_[i2_],
                    "x" => s3_[i2_],
                    "fillcolor" => "#ff1993",
                ),
            ),
        ),
        Dict(
            "showlegend" => false,
            pair_title(s1, "Enrichment = <b>$s2</b>"),
            "yaxis3" =>
                Dict("domain" => (0.328, 1), pair_title("Δ Enrichment")),
            "yaxis2" => Dict(
                "domain" => (0.248, 0.32),
                pair_title("Set"),
                "tickvals" => (),
            ),
            "yaxis" => Dict("domain" => (0, 0.24), pair_title("Score")),
            "xaxis" => Dict(
                pair_title("Feature ($um)"),
                "showspikes" => true,
                "spikemode" => "across",
                "spikedash" => "solid",
                "spikethickness" => -1,
                "spikecolor" => "#27221f",
            ),
        ),
    )

end

function number_enrichment(al, s1_, n1_, st__; u1 = 1, u2 = 1000, pr = 0)

    s2_, n2_ = make_score(s1_, n1_)

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

function read_set(pa)

    di::Dict{String, Vector{String}} = parsefile(pa)

    st_ = collect(keys(di))

    in_ = sortperm(st_)

    st_[in_], collect(values(di))[in_]

end

function make_algorithm(st)

    if st == "S0"

        S0()

    elseif st == "S0a"

        S0a()

    elseif st == "D2"

        D2()

    elseif st == "D2w"

        D2w()

    elseif st == "DD"

        DD()

    end

end

########################################

"""
Run data-rank (single-sample) GSEA.

# Arguments

  - `directory`: Output.
  - `tsv`: Feature-by-sample.
  - `json`: Set-to-features.

# Options

  - `--algorithm`: S0 | S0a | D2 | D2w | DD.
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum set fraction in the data.
  - `--number-of-plots`: The number of top sets to plot.
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

    _, s1_, s2_, N1 = table_part(read_table(tsv))

    s3_, s1__ = read_set(json)

    al = make_algorithm(algorithm)

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

    write_table("$pa.tsv", table_part("Set", s3_, s2_, N2))

    bo_ = map(nu_ -> any(isfinite, nu_), eachrow(N2))

    s4_ = s3_[bo_]

    s2__ = s1__[bo_]

    N3 = N2[bo_, :]

    for in_ in CartesianIndices(N3)[index_extreme(vec(N3), number_of_plots)]

        i1, i2 = Tuple(in_)

        s1 = s4_[i1]

        s2 = s2_[i2]

        s3 = text_2(N3[in_])

        write_enrichment(
            "$pa.$s1.$s2.$s3.html",
            al,
            s1_,
            N1[:, i2],
            s2__[i1],
            Dict(pair_title(s1)),
        )

    end

end

########################################

function table_random(u1, nu, al, s1_, nu_, st__; ke_...)

    N = Matrix{Float64}(undef, length(st__), u1)

    um_ = map(s2_ -> length(intersect(s1_, s2_)), st__)

    seed!(nu)

    @showprogress for nd in 1:u1

        N[:, nd] = number_enrichment(
            al,
            s1_,
            nu_,
            map(u2 -> sample(s1_, u2; replace = false), um_);
            ke_...,
        )

    end

    N

end

function table_random!(um, nu, al, st_, fu, in_, N1, st__; ke_...)

    N2 = Matrix{Float64}(undef, length(st__), um)

    seed!(nu)

    @showprogress for nd in 1:um

        N2[:, nd] = number_enrichment(
            al,
            st_,
            map(nu_ -> make_function(shuffle!(in_), nu_, fu), eachrow(N1)),
            st__;
            ke_...,
        )

    end

    N2

end

########################################

function write_table(di, al, s1_, n1_, s2_, s1__, n2_, N1, u1, s3_)

    u2, u3 = size(N1)

    N2 = Matrix{Float64}(undef, u2, 4)

    N2[:, 1] = n2_

    for i1 in 1:u2

        n6_ = N1[i1, :]

        n1 = mean(filter(<(0), n6_))

        n2 = mean(filter(!<(0), n6_))

        n3 = -n1

        n2_[i1] /= if n2_[i1] < 0

            n3

        else

            n2

        end

        for i2 in 1:u3

            N1[i1, i2] /= if N1[i1, i2] < 0

                n3

            else

                n2

            end

        end

    end

    N2[:, 2] = n2_

    in_, n3_, n4_ = number_significance(n2_, N1)

    N2[in_, 3] = n3_

    N2[in_, 4] = n4_

    write_table(
        joinpath(di, "enrichment.tsv"),
        table_part(
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
        vcat(index_extreme(n5_, u1), filter!(!isnothing, indexin(s3_, s4_))),
    )

        s1 = s4_[nd]

        s2 = text_2(n5_[nd])

        write_enrichment(
            joinpath(di, "$s1.$s2.html"),
            al,
            s1_,
            n1_,
            s2__[nd],
            Dict(pair_title(s1)),
        )

    end

end

########################################

"""
Run user-rank (pre-rank) GSEA.

# Arguments

  - `directory`: Output.
  - `tsv`: Feature-by-metric; the first two columns are features and scores.
  - `json`: Set-to-features.

# Options

  - `--algorithm`: S0 | S0a | D2 | D2w | DD.
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum set fraction in the data.
  - `--number-of-permutations`: For calculating significance.
  - `--seed`: Random seed.
  - `--number-of-plots`: The number of top sets to plot.
  - `--more-plots`: ;-separated set names to plot.
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

    al = make_algorithm(algorithm)

    s1_, nu_ = eachcol(read_table(tsv; select = 1:2))

    s2_, st__ = read_set(json)

    ke_ = (u1 = minimum, u2 = maximum, pr = fraction)

    write_table(
        directory,
        al,
        s1_,
        nu_,
        s2_,
        st__,
        number_enrichment(al, s1_, nu_, st__; ke_...),
        table_random(number_of_permutations, seed, al, s1_, nu_, st__; ke_...),
        number_of_plots,
        split(more_plots, ';'),
    )

end

########################################

"""
Run metric-rank (standard) GSEA.

# Arguments

  - `directory`: Output.
  - `tsv1`: Feature-by-sample; the first two rows are samples and phenotypes.
  - `tsv2`: Feature-by-sample.
  - `json`: Set-to-features.

# Options

  - `--metric`: "signal-to-noise-ratio" | "mean-difference" | "log-ratio".
  - `--algorithm`: S0 | S0a | D2 | D2w | DD.
  - `--minimum`: The minimum set size.
  - `--maximum`: The maximum set size.
  - `--fraction`: The minimum set fraction in the data.
  - `--permutation`: sample | set.
  - `--number-of-permutations`: For calculating significance.
  - `--seed`: Random seed.
  - `--number-of-plots`: The number of top sets to plot.
  - `--more-plots`: ;-separated set names to plot.
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

    _, _, s1_, N1 = table_part(read_table(tsv1))

    in_ = N1[1, :]

    st, s2_, s3_, N1 = table_part(read_table(tsv2))

    N2 = N1[:, indexin(s1_, s3_)]

    fu = if metric == "mean-difference"

        number_difference

    elseif metric == "log-ratio"

        number_ratio

    elseif metric == "signal-to-noise-ratio"

        number_signal

    end

    n1_ = map(
        n2_ -> fu(n2_[findall(isone, in_)], n2_[findall(==(2), in_)]),
        eachrow(N2),
    )

    write_table(
        joinpath(directory, "metric.tsv"),
        table_part(st, s2_, [metric], reshape(n1_, :, 1)),
    )

    al = make_algorithm(algorithm)

    s4_, st__ = read_set(json)

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

            table_random(
                number_of_permutations,
                seed,
                al,
                s2_,
                n1_,
                st__;
                ke_...,
            )

        elseif permutation == "sample"

            table_random!(
                number_of_permutations,
                seed,
                al,
                s2_,
                fu,
                in_,
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
