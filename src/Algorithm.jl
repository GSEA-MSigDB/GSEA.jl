module Algorithm

struct KS0 end

struct A0 end

struct DA2 end

struct DA2W end

struct DA2W0W end

function make(al)

    if al == "ks0"

        KS0()

    elseif al == "a0"

        A0()

    elseif al == "da2"

        DA2()

    elseif al == "da2w"

        DA2W()

    elseif al == "da2w0w"

        DA2W0W()

    end

end

end
