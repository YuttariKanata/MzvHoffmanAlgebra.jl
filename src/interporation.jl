#[ interporation.jl ]#

# This file defines t-MZV calculation.

"""
###################################################################################################
                                        t-MZVs
###################################################################################################
"""

"""
    shift(w::Hoffman; t=1::Rational{BigInt})
    shift(w::Index; t=1::Rational{BigInt})

Computes shift operator S^t(w) of one element.
The shift operator is defined on the admissible Hoffman algebra (words of 0 and 1).
For `Index` elements, they are converted to `Hoffman` elements based on the `orientation`, multiplied, and converted back.
"""

function monomial_shift(w::HoffmanWord; t=1::Rational{BigInt})

end
