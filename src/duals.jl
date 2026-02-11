#[ duals.jl ]#

# This file defines duals, automorphisms, and homomorphisms.

"""
###################################################################################################
                                        Homomorphisms
###################################################################################################
"""

"""
    Hoffman_hom(w::HoffmanWord, image::Vector{Hoffman}) -> Hoffman

Computes the image of a Hoffman word `w` under a homomorphism determined by `image`.
`image[1]` corresponds to the image of `x` (0), and `image[2]` to `y` (1).
"""
function Hoffman_hom(w::HoffmanWord, image::Vector{Hoffman})::Hoffman
    s = one(Hoffman)
    for val in w.t
        # val is 0 or 1. image is 1-based.
        s = s * image[val+1]
    end
    return s
end

"""
    Hoffman_antihom(w::HoffmanWord, image::Vector{Hoffman}) -> Hoffman

Computes the image of a Hoffman word `w` under an anti-homomorphism determined by `image`.
Reverses the order of multiplication.
"""
function Hoffman_antihom(w::HoffmanWord, image::Vector{Hoffman})::Hoffman
    s = one(Hoffman)
    for val in w.t
        s = image[val+1] * s
    end
    return s
end

"""
###################################################################################################
                                        Duals and Automorphisms
###################################################################################################
"""

# Standard Dual (Anti-automorphism)
# Reverses word and swaps x(0) <-> y(1)
function monomial_dual(w::HoffmanWord; kw...)::HoffmanWord
    return HoffmanWord(Tuple(1 - x for x in reverse(w.t)))
end

# Optimized dual for IndexWord (Direct implementation without converting to Hoffman)
# Based on the logic that dual of (k1, ..., kr) corresponds to a specific transformation of partial sums.
function monomial_dual(w::IndexWord; orientation::Symbol=:left)::IndexWord
    if isempty(w) return w end
    
    # Logic for Left orientation (y x^{k-1})
    # This algorithm constructs the dual index directly.
    # Corresponds to monomial_dual_i in oldhoffman.jl
    if orientation === :left
        v = collect(w.t)
        lw = length(v)
        total_weight = sum(v)
        l_dual = total_weight - lw
        
        dual_v = ones(Int, l_dual)
        cid = 0
        for i in lw:-1:1
            cid += v[i] - 1
            if cid > 0 && cid <= l_dual
                dual_v[cid] += 1
            end
        end
        return IndexWord(Tuple(dual_v))
    elseif orientation === :right
        # For right orientation, dual(w) = reverse(dual_left(reverse(w)))
        # Or we can derive a direct algorithm, but reversing is cheap enough for now(2026/2/11).
        return reverse(monomial_dual(reverse(w); orientation=:left))
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end
end

"""
    dual(h::Hoffman) -> Hoffman
    dual(i::Index; orientation::Symbol=:left) -> Index

Computes the standard dual of a Hoffman element or Index element.
For Hoffman words, it reverses the word and swaps x and y.
"""
function dual(h::Hoffman; ks...)::Hoffman
    new_terms = Dict{HoffmanWord, Rational{BigInt}}()
    for (w, c) in h.terms
        wd = monomial_dual(w)
        _add!(new_terms, wd, c)
    end
    return Hoffman(new_terms)
end

function dual(idx::Index; orientation::Symbol=:left)::Index
    new_terms = Dict{IndexWord, Rational{BigInt}}()
    for (w, c) in idx.terms
        _add!(new_terms, monomial_dual(w, orientation=orientation), c)
    end
    return Index(new_terms)
end

# Hoffman Dual (Automorphism)
# Preserves y at start (left) or end (right), swaps rest x <-> y.
function monomial_hof_dual(w::HoffmanWord; orientation::Symbol=:left)::HoffmanWord
    if isempty(w) return w end
    t = w.t
    if orientation === :left
        # Forces first letter to be y (1), swaps the rest
        new_t = (1, (1-x for x in t[2:end])...)
        return HoffmanWord(new_t)
    elseif orientation === :right
        # Forces last letter to be y (1), swaps the rest
        new_t = ((1-x for x in t[1:end-1])..., 1)
        return HoffmanWord(new_t)
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end
end
function monomial_hof_dual(w::IndexWord; orientation::Symbol=:left)::IndexWord
    return IndexWord(monomial_hof_dual(HoffmanWord(w); orientation=orientation)) 
end

"""
    Hoffman_dual(h::Hoffman; orientation::Symbol=:left) -> Hoffman

Computes the Hoffman dual of `h`.
This is an automorphism that fixes `y` (at the "divergent" end) and swaps `x` and `y` elsewhere.
"""
function Hoffman_dual(h::Hoffman; orientation::Symbol=:left)::Hoffman
    new_terms = Dict{HoffmanWord, Rational{BigInt}}()
    for (w, c) in h.terms
        wd = monomial_hof_dual(w, orientation=orientation)
        _add!(new_terms, wd, c)
    end
    return Hoffman(new_terms)
end

function Hoffman_dual(idx::Index; orientation::Symbol=:left)::Index
    # The logic for Hoffman dual on indices is not orientation-independent.
    # The safest and most correct implementation is to convert to Hoffman algebra,
    # apply the well-defined dual there, and convert back.
    return Index(Hoffman_dual(Hoffman(idx; orientation=orientation); orientation=orientation); orientation=orientation)
end

# Landen Dual
function monomial_Landen(w::HoffmanWord; orientation::Symbol=:left)::Hoffman
    s = one(Hoffman)
    hx = Hoffman(x)
    hy = Hoffman(y)
    c1 = hy + hx # x+y
    
    t = w.t
    if orientation === :left
        for val in t
            if val == 1 # y
                s = s * c1
            else # x
                s = -(s * hy)
            end
        end
    elseif orientation === :right
        for i in length(t):-1:1
            val = t[i]
            if val == 1
                s = c1 * s
            else
                s = -(hy * s)
            end
        end
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end
    return s
end

"""
    Landen_dual(h::Hoffman; orientation::Symbol=:left) -> Hoffman

Computes the Landen dual of `h`.
"""
function Landen_dual(h::Hoffman; orientation::Symbol=:left)::Hoffman
    s = zero(Hoffman)
    for (w, c) in h.terms
        term = monomial_Landen(w, orientation=orientation)
        s = s + term * c
    end
    return s
end
Landen_dual(idx::Index; orientation::Symbol=:left)::Index = Index(Landen_dual(Hoffman(idx, orientation=orientation); orientation=orientation))

# Star-word to word
function monomial_star_map(w::HoffmanWord; orientation::Symbol=:left)::Hoffman
    if isempty(w) return one(Hoffman) end
    t = w.t
    
    hx = Hoffman(x)
    hy = Hoffman(y)
    hxy = hx + hy
    
    if orientation === :left
        s = one(Hoffman)
        for val in t[2:end]
            if val == 0 # x
                s = s * hx
            else # y
                s = s * hxy
            end
        end
        return hy * s
    elseif orientation === :right
        s = one(Hoffman)
        for val in t[1:end-1]
            if val == 0
                s = hx * s
            else
                s = hxy * s
            end
        end
        return s * hy
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end
end

"""
    starword_to_word(h::Hoffman; orientation::Symbol=:left) -> Hoffman
    starword_to_word(i::Index; orientation::Symbol=:left) -> Index

Converts a "star-word" (representing a Multiple Zeta Star Value) to a linear combination of standard words.
"""
function starword_to_word(h::Hoffman; orientation::Symbol=:left)::Hoffman
    s = zero(Hoffman)
    for (w, c) in h.terms
        add_term = monomial_star_map(w, orientation=orientation)
        s = s + add_term * c
    end
    return s
end

function starword_to_word(i::Index; orientation::Symbol=:left)::Index
    h = Hoffman(i, orientation=orientation)
    h_res = starword_to_word(h, orientation=orientation)
    return Index(h_res, orientation=orientation)
end