module Information

function make(p1, p2)

    p1 * log2(p1 / p2)

end

function make_thermodynamic_depth(p1, p2)

    make(p1, p2) - make(p2, p1)

end

function make_thermodynamic_breadth(p1, p2)

    make(p1, p2) + make(p2, p1)

end

function make_antisymmetric_kullback_leibler_divergence(p1, p2, p3, p4 = p3)

    make(p1, p3) - make(p2, p4)

end

function make_symmetric_kullback_leibler_divergence(p1, p2, p3, p4 = p3)

    make(p1, p3) + make(p2, p4)

end

function make_jensen_shannon_divergence(p1, p2)

    make_symmetric_kullback_leibler_divergence(p1, p2, (p1 + p2) * 0.5)

end

end
