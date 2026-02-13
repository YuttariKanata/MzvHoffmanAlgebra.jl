#[ arithmetic.jl ]#

# This file defines arithmetic functions

import Base: +, -, *, ^, //, inv

#=
export shift_degree
=#

"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## Arithmetic Functions ###############################################################

# isdefined table
#=

  Addition and Subtraction
--------------------------------------------------------------------------------------------------------------------------
|  Add/Sub       |  NN  |  HoffmanWord  |  Hoffman  |  IndexWord  |  Index  |  Poly NN   |  Poly Hoffman  |  Poly Index  |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  NN            |  NN  |  Hof          |  Hof      |  Idx        |  Idx    |  Poly NN   |  Poly Hof      |  Poly Idx    |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  HofW          |      |  Hof          |  Hof      |  X          |  X      |  Poly Hof  |  Poly Hof      |  X           |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  HoffmanWord   |      |               |  Hof      |  X          |  X      |  Poly Hof  |  Poly Hof      |  X           |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  IndexWord     |      |               |           |  Idx        |  Idx    |  Poly Idx  |  X             |  Poly Idx    |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  Index         |      |               |           |             |  Idx    |  Poly Idx  |  X             |  Poly Idx    |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly NN       |      |               |           |             |         |  Poly NN   |  Poly Hof      |  Poly Idx    |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly Hoffman  |      |               |           |             |         |            |  Poly Hof      |  X           |
|----------------+------+---------------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly Index    |      |               |           |             |         |            |                |  Poly Idx    |
--------------------------------------------------------------------------------------------------------------------------

=#



############################## ADDITIVE INVERSE ##############################

# HoffmanWord
(-)(a::HoffmanWord)::Hoffman = Hoffman(a, -1)
(+)(a::HoffmanWord)::HoffmanWord = a

# Hoffman
function (-)(a::Hoffman)::Hoffman
    return Hoffman(Dict{HoffmanWord, Rational{BigInt}}(w => -c for (w, c) in a.terms))
end
(+)(a::Hoffman)::Hoffman = a

# IndexWord
(-)(a::IndexWord)::Index = Index(a, -1)
(+)(a::IndexWord)::IndexWord = a

# Index
function (-)(a::Index)::Index
    return Index(Dict{IndexWord, Rational{BigInt}}(w => -c for (w, c) in a.terms))
end
(+)(a::Index)::Index = a

# Poly
function (-)(a::Poly{A})::Poly{A} where A
    return Poly{A}(Dict{Int, A}(d => -h for (d, h) in a.terms))
end
function (+)(a::Poly{A})::Poly{A} where A
    return a
end


############################## ADD and SUBTRACT ##############################

# Use promotion to handle different types
# Specific Word arithmetic to avoid recursion in promote
+(a::HoffmanWord, b::HoffmanWord) = Hoffman(a) + Hoffman(b)
-(a::HoffmanWord, b::HoffmanWord) = Hoffman(a) - Hoffman(b)
+(a::IndexWord, b::IndexWord) = Index(a) + Index(b)
-(a::IndexWord, b::IndexWord) = Index(a) - Index(b)

# Mixed types with promote (excluding NN+NN which is Base)
for op in [:(+), :(-)]
    @eval begin
        $op(a::Hoffman, b::HoffmanWord) = $op(promote(a,b)...)
        $op(a::HoffmanWord, b::Hoffman) = $op(promote(a,b)...)
        $op(a::Index, b::IndexWord) = $op(promote(a,b)...)
        $op(a::IndexWord, b::Index) = $op(promote(a,b)...)
        $op(a::NN, b::Union{Hoffman, HoffmanWord}) = $op(promote(a,b)...)
        $op(a::Union{Hoffman, HoffmanWord}, b::NN) = $op(promote(a,b)...)
        $op(a::NN, b::Union{Index, IndexWord}) = $op(promote(a,b)...)
        $op(a::Union{Index, IndexWord}, b::NN) = $op(promote(a,b)...)
    end
end

# Unified Hoffman/Index Addition
function +(a::T, b::T)::T where T <: Union{Hoffman, Index}
    # Optimization: start with the larger dictionary and merge the smaller one
    if length(a.terms) < length(b.terms)
        a, b = b, a
    end
    return T(_add_terms(a.terms, b.terms, +))
end

function -(a::T, b::T)::T where T <: Union{Hoffman, Index}
    return T(_add_terms(a.terms, b.terms, -))
end

########## Poly ##########
+(a::Poly, b::Poly) = +(promote(a,b)...)
-(a::Poly, b::Poly) = -(promote(a,b)...)
+(a::Poly, b) = +(promote(a,b)...)
-(a::Poly, b) = -(promote(a,b)...)
+(a, b::Poly) = +(promote(a,b)...)
-(a, b::Poly) = -(promote(a,b)...)

function +(a::Poly{A}, b::Poly{A})::Poly{A} where A
    if length(a.terms) < length(b.terms)
        a, b = b, a
    end
    return Poly{A}(_add_terms(a.terms, b.terms, +))
end
function -(a::Poly{A}, b::Poly{A})::Poly{A} where A
    return Poly{A}(_add_terms(a.terms, b.terms, -))
end

############################## MULTIPLICATION ##############################

*(a::Hoffman, b::HoffmanWord) = *(promote(a,b)...)
*(a::HoffmanWord, b::Hoffman) = *(promote(a,b)...)
*(a::Union{Hoffman, HoffmanWord}, b::NN) = *(promote(a,b)...)
*(a::NN, b::Union{Hoffman, HoffmanWord}) = *(promote(a,b)...)
*(a::Index, b::IndexWord) = *(promote(a,b)...)
*(a::IndexWord, b::Index) = *(promote(a,b)...)
*(a::Union{Index, IndexWord}, b::NN) = *(promote(a,b)...)
*(a::NN, b::Union{Index, IndexWord}) = *(promote(a,b)...)

# specific
@inline *(a::HoffmanWord, b::HoffmanWord)::HoffmanWord = vcat(a, b)
@inline *(a::IndexWord, b::IndexWord)::IndexWord = vcat(a, b)

# base
function *(a::T, b::T)::T where T <: Union{Hoffman, Index}
    new_terms = Dict{keytype(a.terms), valtype(a.terms)}()
    isempty(a.terms) || isempty(b.terms) && return T(new_terms)
    sizehint!(new_terms, length(a.terms) * length(b.terms))
    for (wa, ca) in a.terms
        for (wb, cb) in b.terms
            w = wa*wb
            new_terms[w] = get(new_terms, w, zero(valtype(a.terms))) + ca*cb
        end
    end
    return T(new_terms)
end



########## Poly ##########

*(a::Poly, b::Poly) = *(promote(a,b)...)
*(a::Poly, b) = *(promote(a,b)...)
*(a, b::Poly) = *(promote(a,b)...)

function *(a::Poly{A}, b::Poly{A})::Poly{A} where A

    # 単項最適化
    if is_monomial(a)
        (d,c) = first(a.terms)
        if isone(c)
            return shift_degree(b, d)
        end
    end
    if is_monomial(b)
        (d,c) = first(b.terms)
        if isone(c)
            return shift_degree(a, d)
        end
    end

    # 通常の積
    new_terms = Dict{Int, A}()
    sizehint!(new_terms, length(a.terms) * length(b.terms))
    for (da, ha) in a
        for (db, hb) in b
            d = da + db
            new_terms[d] = get(new_terms, d, zero(A)) + ha*hb
        end
    end
    return Poly{A}(new_terms)
end


############################## POWER ##############################
import Base: literal_pow, power_by_squaring

# HoffmanWord
function ^(t::W, n::Integer)::W where W <: AbstractWord
    if n < 0
        throw(DomainError(n, "negative power is not supported for AbstractWord"))
    elseif n == 0
        return W()
    end
    len = length(t)
    total = len * n
    v = Vector{eltype(t)}(undef, total)
    for i in 0:n-1
        copyto!(v, i*len + 1, t, 1, len)
    end
    return W(v)
end

# Hoffman, Index
function ^(a::Union{Hoffman, Index}, n::Integer)::Union{Hoffman, Index}
    return power_by_squaring(a, n)
end

# Poly
function ^(a::Poly{A}, n::Integer)::Poly{A} where A
    if n < 0
        throw(DomainError(n, "negative power is not supported for Poly"))
    end

    # 最適化: a が単項 T^k (係数が exactly one(Hoffman)) の場合
    if is_monomial(a)
        (deg, coeff) = first(a.terms)
        # Use ^ for coefficient, which will use power_by_squaring for Hoffman/Index
        terms = Dict(deg * Int(n) => coeff^n)
        return Poly{A}(terms)
    end

    return power_by_squaring(a, n)
end

# Literal power optimizations
# x^2, x^3 etc. will be lowered to Base.literal_pow(^, x, Val(2)) etc.
literal_pow(::typeof(^), x::Union{Hoffman, Index, Poly}, ::Val{0}) = one(x)
literal_pow(::typeof(^), x::Union{Hoffman, Index, Poly}, ::Val{1}) = copy(x)
literal_pow(::typeof(^), x::Union{Hoffman, Index, Poly}, ::Val{2}) = x*x
literal_pow(::typeof(^), x::Union{Hoffman, Index, Poly}, ::Val{3}) = x*x*x

# degree shift
""" r , n -> r T^n """
function shift_degree(r::Poly{A},n::Int)::Poly{A} where A
    return Poly{A}(Dict{Int, A}(d+n => h for (d, h) in r.terms))
end

# Inverse (for literal x^-1)
inv(a::Hoffman) = throw(ArgumentError("Hoffman-type powers for negative exponents are not defined"))
inv(a::Index) = throw(ArgumentError("Index-type powers for negative exponents are not defined"))
inv(a::Poly{A}) where A = throw(DomainError(-1, "negative power is not supported for Poly"))


############################# DIVISION for NN ##############################
# NN Word
function //(b::HoffmanWord, a::NN)::Hoffman
    return Hoffman(Dict{HoffmanWord, Rational{BigInt}}(b => 1//a))
end

# NN Hoffman
function //(b::Hoffman, a::NN)::Hoffman
    return Hoffman(Dict{HoffmanWord, Rational{BigInt}}(w => c//a for (w,c) in b.terms))
end

# NN IndexWord
function //(b::IndexWord, a::NN)::Index
    return Index(Dict{IndexWord, Rational{BigInt}}(b => 1//a))
end

# NN Index
function //(b::Index, a::NN)::Index
    return Index(Dict{IndexWord, Rational{BigInt}}(w => c//a for (w,c) in b.terms))
end

# NN Poly
function //(b::Poly{A}, a::NN)::Poly{A} where A
    return Poly{A}(Dict{Int, A}(d => h//a for (d,h) in b.terms))
end
