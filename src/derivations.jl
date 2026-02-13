#[ derivations.jl ]#

# This file defines derivations on the algebras.

"""
###################################################################################################
                                        Derivations
###################################################################################################
"""

"""
    Hoffman_derivation(w::HoffmanWord, image::Vector{Hoffman}) -> Hoffman

Applies a derivation to a word `w`. The derivation is defined by its action on the generators `x` and `y`,
given by `image[1]` and `image[2]` respectively.
"""
function Hoffman_derivation(w::HoffmanWord, image::Vector{Hoffman})::Hoffman
    s = zero(Hoffman)
    lw = length(w)
    for i in 1:lw
        # Leibniz rule: d(abc) = d(a)bc + ad(b)c + abd(c)
        # Here, we apply the derivation to the i-th letter.
        prefix = Hoffman(w[1:i-1])
        suffix = Hoffman(w[i+1:lw])
        # w[i] is 0 or 1, image is 1-based.
        d_letter = image[w[i]+1]
        
        s = s + prefix * d_letter * suffix
    end
    return s
end

"""
    dell(h::Hoffman, n::Int; orientation::Symbol=:left) -> Hoffman

Computes the derivation `∂_n` on a Hoffman element `h`.
"""
function dell(h::Hoffman, n::Int; orientation::Symbol=:left)::Hoffman
    s = zero(Hoffman)
    # The derivation ∂_n is defined by its action on the generators x and y.
    # This implementation follows the definition from oldhoffman.jl.
    local xyn1::Hoffman
    if orientation === :left
        # ∂_n(x) = y * (x+y)^(n-1) * x
        xyn1 = y * (x+y)^(n-1) * x
    elseif orientation === :right
        # ∂_n(x) = x * (x+y)^(n-1) * y
        xyn1 = x * (x+y)^(n-1) * y
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end
    image = [xyn1, -xyn1]
    for (w, c) in h.terms
        s = s + Hoffman_derivation(w, image) * c
    end
    return s
end

dell(idx::Index, n::Int; orientation::Symbol=:left)::Index = Index(dell(Hoffman(idx; orientation=orientation), n; orientation=orientation); orientation=orientation)
dell(w::HoffmanWord, n::Int; kw...)::Hoffman = dell(Hoffman(w), n; kw...)
dell(w::IndexWord, n::Int; kw...)::Index = Index(dell(Hoffman(w), n; kw...))
