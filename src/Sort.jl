module Sort

function make(an_, nu_)

    id_ = findall(!isnan, nu_)

    an_ = an_[id_]

    nu_ = nu_[id_]

    sortperm!(id_, nu_; rev = true)

    an_[id_], nu_[id_]

end

end
