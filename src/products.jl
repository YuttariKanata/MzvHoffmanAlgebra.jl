#[ products.jl ]#

# This file defines algebraic product operations.

"""
###################################################################################################
                                        Algebraic Operations
###################################################################################################
"""

"""
    shuffle_product(a::Hoffman, b::Hoffman)
    shuffle_product(a::Index, b::Index; orientation::Symbol=:left)

Computes the shuffle product of two elements.
The shuffle product is defined on the Hoffman algebra (words of 0 and 1).
For `Index` elements, they are converted to `Hoffman` elements based on the `orientation`, multiplied, and converted back.

Aliases: `ш`, `⨝`
"""
function shuffle_product(a::Hoffman, b::Hoffman; kw...)::Hoffman
    new_terms = Dict{HoffmanWord, Rational{BigInt}}()
    for (wa, ca) in a.terms
        va = collect(wa)
        for (wb, cb) in b.terms
            vb = collect(wb)
            coeff = ca * cb
            for v in monomial_sh(va, vb)
                _add!(new_terms, HoffmanWord(v), coeff)
            end
        end
    end
    return Hoffman(new_terms)
end


"""
    shuffle_product(a::Index, b::Index; orientation::Symbol=:left)

Optimized shuffle product for Index using monomial_sh_i.
"""
function shuffle_product(a::Index, b::Index; kw...)::Index
    # Dispatch to fast implementation
    new_terms = Dict{IndexWord, Rational{BigInt}}()
    
    for (wa, ca) in a.terms
        va = collect(wa)
        for (wb, cb) in b.terms
            vb = collect(wb)
            coeff = ca * cb
            for (w_res, c_res) in monomial_sh_i(va, vb)
                _add!(new_terms, w_res, coeff * c_res)
            end
        end
    end
    return Index(new_terms)
end

"""
    stuffle_product(a::Index, b::Index)
    stuffle_product(a::Hoffman, b::Hoffman; orientation::Symbol=:left)

Computes the stuffle (harmonic) product of two elements.
The stuffle product is defined on the Index algebra (words of integer indices).
For `Hoffman` elements, they are converted to `Index` elements based on the `orientation`, multiplied, and converted back.

Alias: `∗`
"""
function stuffle_product(a::Index, b::Index; kw...)::Index
    new_terms = Dict{IndexWord, Rational{BigInt}}()
    for (wa, ca) in a.terms
        va = collect(wa)
        for (wb, cb) in b.terms
            vb = collect(wb)
            coeff = ca * cb
            for v in monomial_st(va, vb)
                _add!(new_terms, IndexWord(v), coeff)
            end
        end
    end
    return Index(new_terms)
end

function stuffle_product(a::Hoffman, b::Hoffman; orientation::Symbol=:left)::Hoffman
    ia = Index(a, orientation=orientation)
    ib = Index(b, orientation=orientation)
    i_res = stuffle_product(ia, ib)
    return Hoffman(i_res, orientation=orientation)
end

"""
    star_stuffle_product(a::Index, b::Index; orientation::Symbol=:left)
    star_stuffle_product(a::Hoffman, b::Hoffman; orientation::Symbol=:left)

Computes the star-stuffle product of two elements.
This corresponds to the product of Multiple Zeta Star Values (MZSV).

Alias: `⋆`
"""
function star_stuffle_product(a::Index, b::Index; kw...)::Index
    new_terms = Dict{IndexWord, Rational{BigInt}}()
    for (wa, ca) in a.terms
        va = collect(wa)
        la = length(va)
        for (wb, cb) in b.terms
            vb = collect(wb)
            lb = length(vb)
            
            parity_total = (la + lb) & 1
            coeff = ca * cb
            
            for v in monomial_st(va, vb)
                # Sign is (-1)^(dep(u) + dep(v) - dep(w))
                # If parity matches, sign is +, else -
                if (length(v) & 1) == parity_total
                    _add!(new_terms, IndexWord(v), coeff)
                else
                    _add!(new_terms, IndexWord(v), -coeff)
                end
            end
        end
    end
    return Index(new_terms)
end

function star_stuffle_product(a::Hoffman, b::Hoffman; orientation::Symbol=:left)::Hoffman
    ia = Index(a, orientation=orientation)
    ib = Index(b, orientation=orientation)
    i_res = star_stuffle_product(ia, ib)
    return Hoffman(i_res, orientation=orientation)
end

for op in [:(shuffle_product), :(stuffle_product), :(star_stuffle_product)]
    @eval begin
        $op(a::Hoffman, b::HoffmanWord; kw...)     = $op(Hoffman(a), Hoffman(b); kw...)
        $op(a::HoffmanWord, b::Hoffman; kw...)     = $op(Hoffman(a), Hoffman(b); kw...)
        $op(a::HoffmanWord, b::HoffmanWord; kw...) = $op(Hoffman(a), Hoffman(b); kw...)
        $op(a::Index, b::IndexWord; kw...)         = $op(Index(a), Index(b); kw...)
        $op(a::IndexWord, b::Index; kw...)         = $op(Index(a), Index(b); kw...)
        $op(a::IndexWord, b::IndexWord; kw...)     = $op(Index(a), Index(b); kw...)
    end
end

# Squaring optimizations (Internal use)
# Computes a * a efficiently by exploiting commutativity/symmetry where possible.
function shuffle_square(a::Hoffman)::Hoffman
    new_terms = Dict{HoffmanWord, Rational{BigInt}}()
    terms = collect(a.terms)
    n = length(terms)
    for i in 1:n
        (wi, ci) = terms[i]
        vi = collect(wi)
        # Diagonal term: ci^2 * (wi ш wi)
        # Note: monomial_sh(vi, vi) is symmetric, but we use the standard function.
        coeff_sq = ci^2
        for v in monomial_sh(vi, vi)
            _add!(new_terms, HoffmanWord(v), coeff_sq)
        end
        
        # Cross terms: 2 * ci * cj * (wi ш wj)
        for j in i+1:n
            (wj, cj) = terms[j]
            vj = collect(wj)
            coeff_cross = 2 * ci * cj
            for v in monomial_sh(vi, vj)
                _add!(new_terms, HoffmanWord(v), coeff_cross)
            end
        end
    end
    return Hoffman(new_terms)
end

function shuffle_square(a::Index; kw...)::Index
    return shuffle_product(a,a)
end

function stuffle_square(a::Index)::Index
    new_terms = Dict{IndexWord, Rational{BigInt}}()
    terms = collect(a.terms)
    n = length(terms)
    for i in 1:n
        (wi, ci) = terms[i]
        vi = collect(wi)
        coeff_sq = ci^2
        for v in monomial_st(vi, vi)
            _add!(new_terms, IndexWord(v), coeff_sq)
        end
        
        for j in i+1:n
            (wj, cj) = terms[j]
            vj = collect(wj)
            coeff_cross = 2 * ci * cj
            for v in monomial_st(vi, vj)
                _add!(new_terms, IndexWord(v), coeff_cross)
            end
        end
    end
    return Index(new_terms)
end

function stuffle_square(a::Hoffman; orientation::Symbol=:left)::Hoffman
    ia = Index(a, orientation=orientation)
    i_res = stuffle_square(ia)
    return Hoffman(i_res, orientation=orientation)
end

function star_stuffle_square(a::Index; kw...)::Index
    # Similar logic to star_stuffle_product but optimized for square
    # Since we don't have a specialized monomial_st_star, we use the product logic with square optimization loops.
    new_terms = Dict{IndexWord, Rational{BigInt}}()
    terms = collect(a.terms)
    n = length(terms)
    
    # Helper to compute star stuffle term for two words
    function add_star_term!(d, va, vb, coeff)
        words = monomial_st(va, vb)
        parity_total = (length(va) + length(vb)) & 1
        for v in words
            if (length(v) & 1) == parity_total
                _add!(d, IndexWord(v), coeff)
            else
                _add!(d, IndexWord(v), -coeff)
            end
        end
    end

    for i in 1:n
        (wi, ci) = terms[i]
        vi = collect(wi)
        add_star_term!(new_terms, vi, vi, ci^2)
        
        for j in i+1:n
            (wj, cj) = terms[j]
            vj = collect(wj)
            add_star_term!(new_terms, vi, vj, 2 * ci * cj)
        end
    end
    return Index(new_terms)
end

function star_stuffle_square(a::Hoffman; orientation::Symbol=:left)::Hoffman
    ia = Index(a, orientation=orientation)
    i_res = star_stuffle_square(ia)
    return Hoffman(i_res, orientation=orientation)
end

for op in [:(shuffle_square), :(stuffle_square), :(star_stuffle_square)]
    @eval begin
        $op(a::HoffmanWord; kw...) = $op(Hoffman(a); kw...)
        $op(a::IndexWord; kw...) = $op(Index(a); kw...)
    end
end

# Power functions (Binary Exponentiation)
# Generic helper for binary exponentiation
function _bin_pow(a, n::Integer, op, sq_op)
    if n < 0
        throw(DomainError(n, "Negative power not supported"))
    elseif n == 0
        return one(typeof(a))
    elseif n == 1
        return a
    end
    
    base = a
    
    # Handle trailing zeros (initial squarings)
    t = trailing_zeros(n)
    for _ in 1:t
        base = sq_op(base)
    end
    n >>= t
    
    # First bit is now 1
    res = base
    n >>= 1
    
    while n > 0
        base = sq_op(base)
        if (n & 1) == 1
            res = op(res, base)
        end
        n >>= 1
    end
    return res
end

# Iterative power (naive) - faster for sparse bases
function _iter_pow(a, n::Integer, op)
    if n == 0 return one(typeof(a)) end
    res = a
    for _ in 2:n
        res = op(res, a)
    end
    return res
end

"""
    shuffle_pow(a, n::Integer)

Computes the n-th power of `a` with respect to the shuffle product.
`a ш a ш ... ш a`
"""
function shuffle_pow(a::Hoffman, n::Integer; kw...)
    # Heuristic: if a is sparse (e.g. x+y), iterative multiplication is faster
    if length(a) <= 2
        return _iter_pow(a, n, shuffle_product)
    end
    _bin_pow(a, n, shuffle_product, shuffle_square)
end

function shuffle_pow(a::Index, n::Integer; kw...)

    if length(a) <= 2
        return _iter_pow(a, n, shuffle_product)
    end
    _bin_pow(a, n, shuffle_product, shuffle_square)
end

# Stuffle power (similar logic)
"""
    stuffle_pow(a, n::Integer)

Computes the n-th power of `a` with respect to the stuffle product.
`a ∗ a ∗ ... ∗ a`
"""
function stuffle_pow(a::Index, n::Integer; kw...)
    if length(a) <= 2
        return _iter_pow(a, n, stuffle_product)
    end
    _bin_pow(a, n, stuffle_product, stuffle_square)
end

function stuffle_pow(a::Hoffman, n::Integer; orientation::Symbol=:left)
    op = (x, y) -> stuffle_product(x, y, orientation=orientation)
    sq_op = (x) -> stuffle_square(x, orientation=orientation)
    
    if length(a) <= 2
        return _iter_pow(a, n, op)
    end
    _bin_pow(a, n, op, sq_op)
end

"""
    star_stuffle_pow(a, n::Integer)

Computes the n-th power of `a` with respect to the star-stuffle product.
`a ⋆ a ⋆ ... ⋆ a`
"""
function star_stuffle_pow(a::Index, n::Integer; kw...)
    if length(a) <= 2
        return _iter_pow(a, n, star_stuffle_product)
    end
    _bin_pow(a, n, star_stuffle_product, star_stuffle_square)
end

function star_stuffle_pow(a::Hoffman, n::Integer; orientation::Symbol=:left)
    op = (x, y) -> star_stuffle_product(x, y, orientation=orientation)
    sq_op = (x) -> star_stuffle_square(x, orientation=orientation)
    
    if length(a) <= 2
        return _iter_pow(a, n, op)
    end
    _bin_pow(a, n, op, sq_op)
end

for op in [:(shuffle_pow), :(stuffle_pow), :(star_stuffle_pow)]
    @eval begin
        $op(a::HoffmanWord, n::Integer; kw...) = $op(Hoffman(a), n; kw...)
        $op(a::IndexWord, n::Integer; kw...) = $op(Index(a), n; kw...)
    end
    
end

"""
    st_index1_pow(n::Int) -> Index

Computes the stuffle product of Index(1) with itself n times: (1) * ... * (1).
This is equal to sum_{compositions of n} (n! / prod(parts!)) * Index(composition).
Used for regularization polynomials.
"""
function st_index1_pow(n::Int)::Index
    if n < 0
        throw(DomainError(n, "Negative power not supported"))
    elseif n == 0
        return one(Index)
    end

    idx = zero(Index)
    nf = factorial(BigInt(n))
    
    # Iterate over all compositions of n by iterating over Hoffman words of length n starting with 1.
    limit = 1 << (n-1)
    bits = Vector{Int}(undef, n)
    bits[1] = 1
    
    for mask in 0:(limit-1)
        for i in 0:(n-2)
            bits[i+2] = (mask >> (n-2-i)) & 1
        end
        
        hw = HoffmanWord(bits)
        iv = idxprs(hw)
        iw = IndexWord(iv)
        
        coeff = multinomial(iv, nf)
        _add!(idx.terms, iw, Rational{BigInt}(coeff))
    end
    
    return idx
end

"""
    sh_y_pow(n::Int) -> Hoffman

Computes the shuffle product of y with itself n times: y ш ... ш y.
This is equal to n! * y^n.
Used for regularization polynomials.
"""
function sh_y_pow(n::Int)::Hoffman
    if n < 0
        throw(DomainError(n, "Negative power not supported"))
    elseif n == 0
        return one(Hoffman)
    end
    
    # y^n corresponds to HoffmanWord of all 1s with length n
    w = HoffmanWord(ones(Int, n))
    c = factorial(BigInt(n))
    
    return Hoffman(Dict(w => Rational{BigInt}(c)))
end

const ш = shuffle_product
const ⨝ = shuffle_product
const ∗ = stuffle_product
const ⋆ = star_stuffle_product