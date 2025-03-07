module Algorithm

struct KS end

struct KSA end

struct KLIOM end

struct KLIOP end

struct KLI end

struct KLI1 end

function make(al)

    if al == "ks"

        KS()

    elseif al == "ksa"

        KSA()

    elseif al == "kliom"

        KLIOM()

    elseif al == "kliop"

        KLIOP()

    elseif al == "kli"

        KLI()

    elseif al == "kli1"

        KLI1()

    end

end

end
