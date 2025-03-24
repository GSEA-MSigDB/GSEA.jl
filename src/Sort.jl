module Sort

function make(an_, nu_)

    in_ = findall(!isnan, nu_)

    an_ = an_[in_]

    nu_ = nu_[in_]

    sortperm!(in_, nu_; rev = true)

    an_[in_], nu_[in_]

end

end
