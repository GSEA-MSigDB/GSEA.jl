module Algorithm

struct KS0 end

struct A0 end

struct DA2 end

struct DA2W end

struct DA2W0W end

function make(st)

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

end
