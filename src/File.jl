module File

using Nucleus

function read_cls(cl)

    l1, l2, l3 = readlines(cl)

    l2 = l2[2:end]

    s3_ = split(l3)

    st_, N = if l1 == "#numeric"

        l2, [parse(Float64, st) for _ in 1:1, st in s3_]

    else

        s1_ = split(l1)

        s2_ = split(l2)

        @assert parse(Int, s1_[1]) == lastindex(s3_)

        @assert parse(Int, s1_[2]) == lastindex(s2_) == lastindex(unique(s3_))

        di = Dict(s2_[id] => id for id in eachindex(s2_))

        join(s2_, '_'), [di[st] for _ in 1:1, st in s3_]

    end

    "Phenotype", st_, map!(id -> "Sample $id", s3_, eachindex(s3_)), N

end

function read_gct(gc)

    Nucleus.Table.ge(Nucleus.Table.rea(gc; header = 3, drop = ["Description"]))

end

function read_gmt(gm)

    di = Dict{String, Vector{String}}()

    for li in eachline(gm)

        st_ = split(li, '\t')

        st = st_[1]

        @assert !haskey(di, st)

        di[st] = filter!(!isempty, st_[3:end])

    end

    di

end

end
