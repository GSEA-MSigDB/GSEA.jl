module File

using Nucleus

function read_cls(cl)

    l1, l2, l3 = readlines(cl)

    na = "Phenotype"

    l2 = l2[2:end]

    s3_ = split(l3)

    na_ = map(id -> "Sample $id", eachindex(s3_))

    if l1 == "#numeric"

        Nucleus.Table.make(na, l2, na_, [parse(Float64, st) for _ in 1:1, st in s3_])

    else

        s1_ = split(l1)

        s2_ = split(l2)

        @assert parse(Int, s1_[1]) == lastindex(s3_)

        @assert parse(Int, s1_[2]) == lastindex(s2_) == lastindex(unique(s3_))

        di = Dict(st => id for (id, st) in enumerate(s2_))

        Nucleus.Table.make(na, join(s2_, '_'), na_, [di[st] for _ in 1:1, st in s3_])

    end

end

function read_gct(gc)

    Nucleus.Table.rea(gc; header = 3, drop = ["Description"])

end

function read_gmt(gm)

    di = Dict{String, Vector{String}}()

    for li in eachline(gm)

        sp_ = split(li, '\t')

        na = sp_[1]

        @assert !haskey(di, na)

        di[na] = filter!(!isempty, sp_[3:end])

    end

    di

end

end
