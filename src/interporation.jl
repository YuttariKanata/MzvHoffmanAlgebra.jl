#[ interporation.jl ]#

# This file defines t-MZV calculation.

"""
###################################################################################################
                                        t-MZVs
###################################################################################################
"""

"""
    shift(w::Hoffman; t=0::Rational{BigInt})::Hoffman
    shift(idx::Index; t=0::Rational{BigInt})::Index

Computes shift operator S^t(w) of one element.
The shift operator is defined on the admissible Hoffman algebra (words of 0 and 1).
For `Index` elements, they are converted to `Hoffman` elements based on the `orientation`, multiplied, and converted back.
"""

function monomial_shift(w::HoffmanWord; t=0//1::Rational{BigInt})
    k = HoffmanWord(w.t[2:end])
    return y*Hoffman_hom(k, [x, t*x + y])
end

function monomial_shift(w::IndexWord; t=0//1::Rational{BigInt})
    k = HoffmanWord(HoffmanWord(w).t[2:end])
    return Index(y*Hoffman_hom(k, [x, t*x + y]))
end

function shift(w::Hoffman; t=0//1::Rational{BigInt})
    ans = zero(Hoffman)
    for (h, c) in w.terms
        hd = c*monomial_shift(h; t=t)
        ans += hd
    end
    return ans
end

shift(w::HoffmanWord; t = 0//1::Rational{BigInt}) = monomial_shift(w; t=t)
shift(idx::Index; t=0//1::Rational{BigInt}) = Index(shift(Hoffman(idx); t=t))
shift(idx::IndexWord; t=0//1::Rational{BigInt}) = monomial_shift(idx; t=t)

"""
    sharp_shift(w::Hoffman)
    sharp_shift(idx::Index)

Computes "sharp-MZV" such that 2^weight * shift(;t=1/2).
The shift operator is defined on the admissible Hoffman algebra (words of 0 and 1).
For `Index` elements, they are converted to `Hoffman` elements based on the `orientation`, multiplied, and converted back.

sharp_shifted Indexes are important e.g. Zhao's 2-1 formula.
"""

monomial_sharp_shift(kk::HoffmanWord) = 2^sum(kk.t)*shift(kk; t=big(1)//2)
monomial_sharp_shift(kk::IndexWord) = 2^length(kk.t)*shift(kk; t=big(1)//2)

function sharp_shift(w::Hoffman)
    ans = zero(Hoffman)
    for (h, c) in w.terms
        hd = c*monomial_sharp_shift(h)
        ans += hd
    end
    return ans
end

function sharp_shift(idx::Index)
    ans = zero(Hoffman)
    for (h, c) in idx.terms
        hd = c*monomial_sharp_shift(h)
        ans += hd
    end
    return ans
end

sharp_shift(kk::HoffmanWord) = monomial_sharp_shift(kk)
sharp_shift(kk::IndexWord) = monomial_sharp_shift(kk)
