#[ converting.jl ]#

# This file defines conversion functions between each type

#=
export index, x, y
=#

"""
###################################################################################################
                                        Operators
###################################################################################################
"""


###################################################################################################
############## Conversion Functions ###############################################################


# Operator()に通す -> 正規化される運命である ことの理解
# だから.opsを無理にいじることは推奨されていない(unsafe)

# 正規化される前を扱うコンストラクタ
# Operator() is defined in types.jl

function Operator(op::AbstractOp)::Operator
    # An Operator is immutable, build the vector first.
    if op.cnt == 0
        return Operator()
    end
    # The constructor calls clean(), which also copies.
    return Operator([op])
end

# 最後にOperator()を通すことで正規化するコンストラクタ
function Operator(sym::AbstractString, n::Integer=1)::Operator
    if n == 0
        return Operator()
    end
    return Operator(str_to_op(sym,n))   # 正規化される
end
function Operator(v::Vector{<:Tuple{AbstractString, Integer}})::Operator
    regv = Vector{AbstractOp}(undef,lastindex(v))
    for (i, vi) in enumerate(v)
        regv[i] = str_to_op(vi) # UpDownの正規化がされていない
    end
    return Operator(regv)   # 正規化される
end
Operator(op::Operator)::Operator = op # It's immutable, no copy needed
Operator(op::Type{<:AbstractOp})::Operator = Operator(op(1))
function Operator(op::Type{<:AbstractOp},cnt::Int,n::Int=1)::Operator
    if cnt == 0
        return Operator()
    end
    if op === OpDeriv
        return Operator(OpDeriv(cnt,n))
    else
        return Operator(op(cnt))  # 正規化される
    end
end





"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## Conversion Functions ###############################################################

# [============== about MonoIndex ==============]
IndexWord()::IndexWord = IndexWord(Tuple{}())
IndexWord(v::Vector)::IndexWord = IndexWord(Tuple{Vararg{Int}}(v))
function IndexWord(i::Index)::IndexWord
    if !is_monomial(i)
        throw(DomainError(i,"i must be monomial"))
    end
    return IndexWord(first(keys(i)))
end
function IndexWord(w::Hoffman; orientation::Symbol = :left)::IndexWord
    if !is_monomial(w)
        throw(DomainError(w,"w must be monomial"))
    end
    return IndexWord(first(keys(w)); orientation=orientation)
end
function IndexWord(w::HoffmanWord; orientation::Symbol = :left)::IndexWord
    if orientation === :left # z_k = y*x^{k-1}
        if isempty(w) || w[1] != 1 # yに対応するだけ
            throw(DomainError(w,"w does not start with y"))
        end
        return IndexWord(idxprs(w))
    elseif orientation === :right # z_k = x^{k-1}*y
        if isempty(w) || w[end] != 1
            throw(DomainError(w,"w does not end with y"))
        end
        return IndexWord(idxprs_r(w))
    end
    throw(ArgumentError("orientation must be :left or :right"))
end

IndexWord(a::Int...)::IndexWord = IndexWord(a)
IndexWord(m::IndexWord)::IndexWord = m

index(s...)::IndexWord = IndexWord(s...)    # = Index(s...)にするべき？
# TODO: MonoIndex for hrm shf mpl 


# [============== about Word ==============]
HoffmanWord()::HoffmanWord = HoffmanWord(Tuple{}())
HoffmanWord(s::Int...)::HoffmanWord = HoffmanWord(s)
const x = HoffmanWord(0)
const y = HoffmanWord(1)

HoffmanWord(v::Vector{Int})::HoffmanWord = HoffmanWord(Tuple{Vararg{Int}}(v))
function HoffmanWord(w::Hoffman; orientation::Symbol = :left)::HoffmanWord
    if !is_monomial(w)
        throw(DomainError(w,"w must be monomial"))
    end
    p = first(pairs(w))
    if !isone(p.second)
        throw(DomainError(w,"w's coefficient is not 1"))
    end
    return HoffmanWord(p.first; orientation=orientation)
end

function HoffmanWord(w::IndexWord; orientation::Symbol = :left)::HoffmanWord
    if orientation === :left
        return idxdprs(collect(w))
    elseif orientation === :right
        return idxdprs_r(collect(w))
    end
    throw(ArgumentError("orientation must be :left or :right"))
end
function HoffmanWord(i::Index; orientation::Symbol = :left)::HoffmanWord
    if !is_monomial(i)
        throw(DomainError(i,"i must be monomial"))
    end
    p = first(pairs(i))
    if !isone(p.second)
        throw(DomainError(i,"i's coefficient is not 1"))
    end
    return HoffmanWord(p.first; orientation=orientation)
end

HoffmanWord(w::HoffmanWord)::HoffmanWord = w # Immutable
# TODO: Word for hrm shf mpl


# [============== about Index ==============]
Index()::Index = Index(Dict{IndexWord,Rational{BigInt}}())
function Index(idx::Int...)::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(idx) => Rational{BigInt}(1) ) )
end
function Index(t::Tuple{Vararg{Int}})::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(t) => Rational{BigInt}(1) ) )
end
function Index(v::Vector{Int})::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(v) => Rational{BigInt}(1) ) )
end
function Index(c::NN, v::Vector{Int})::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(v) => Rational{BigInt}(c) ) )
end
function Index(v::Vector{Vector{Int}})::Index
    idx = Index()
    c = Rational{BigInt}(1)
    for vi in v
        wvi = IndexWord(vi)
        idx.terms[wvi] = getindex(idx,wvi) + c  # 1を足すだけなので0になる心配がない
    end
    return idx
end

# star-stuffle用 (Boolは+-がはいる)
function Index(v::Tuple{Vector{Vector{Int}}, Vector{Bool}})::Index
    idx = Index()
    c = Rational{BigInt}(1)
    for i in 1:lastindex(v[1])
        wvi = IndexWord(v[1][i])
        if v[2][i]
            idx.terms[wvi] = getindex(idx,wvi) - c
        else
            idx.terms[wvi] = getindex(idx,wvi) + c
        end
    end
    filter!(p->!iszero(p.second),idx.terms) # 1を足すだけではないので0になる心配がある
    return idx
end
function Index(m::IndexWord,coeff::NN = Rational{BigInt}(1))::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( m => Rational{BigInt}(coeff) ) )
end
function Index(w::HoffmanWord,coeff::NN = Rational{BigInt}(1); orientation::Symbol = :left)::Index
    return Index( Dict{IndexWord,Rational{BigInt}}( IndexWord(w; orientation=orientation) => Rational{BigInt}(coeff) ) )
end
function Index(w::Hoffman; orientation::Symbol = :left)::Index   # Hoffmanは正規化されている
    terms = Dict{IndexWord, Rational{BigInt}}(
        IndexWord(k; orientation=orientation) => c for (k, c) in w.terms
    )
    return Index(terms)
end

Index(i::Index)::Index = i # Immutable
# function Index(n::NN)::Index
#     idx = Index()
#     idx.terms[Word()] = n
#     return idx
# end
# TODO: Index for hrm shf mpl


# [============== about Hoffman ==============]
Hoffman()::Hoffman = Hoffman(Dict{HoffmanWord,Rational{BigInt}}())
function Hoffman(v::Vector{Vector{Int}})::Hoffman
    w = Hoffman()
    c = Rational{BigInt}(1)
    for vi in v
        wvi = HoffmanWord(vi)
        w.terms[wvi] = getindex(w,wvi) + c  # 1を足すだけなので0になる心配がない
    end
    return w
end
function Hoffman(v::Vector{Int})::Hoffman
    return Hoffman( Dict{HoffmanWord,Rational{BigInt}}( HoffmanWord(v) => Rational{BigInt}(1) ) )
end
function Hoffman(m::IndexWord, coeff::NN=Rational{BigInt}(1); orientation::Symbol = :left)::Hoffman
    return Hoffman( Dict{HoffmanWord,Rational{BigInt}}( HoffmanWord(m; orientation=orientation) => Rational{BigInt}(coeff) ) )
end
function Hoffman(wm::HoffmanWord, coeff::NN=Rational{BigInt}(1))::Hoffman
    return Hoffman( Dict{HoffmanWord,Rational{BigInt}}( wm => Rational{BigInt}(coeff) ) )
end
function Hoffman(i::Index; orientation::Symbol = :left)::Hoffman # Indexは正規化されている
    terms = Dict{HoffmanWord, Rational{BigInt}}(
        HoffmanWord(k; orientation=orientation) => c for (k, c) in i.terms
    )
    return Hoffman(terms)
end

Hoffman(w::Hoffman)::Hoffman = w # Immutable
# function Hoffman(n::NN)::Hoffman
#     w = Hoffman()
#     w.terms[Word()] = n
#     return w
# end
# TODO: Hoffman for hrm shf mpl


# [============== about Poly ==============]

import Base: convert, promote_rule

# --- Promotion Rules ---
promote_rule(::Type{Hoffman}, ::Type{<:Union{HoffmanWord, NN}}) = Hoffman
promote_rule(::Type{Index}, ::Type{<:Union{IndexWord, NN}}) = Index
promote_rule(::Type{HoffmanWord}, ::Type{<:NN}) = Hoffman
promote_rule(::Type{IndexWord}, ::Type{<:NN}) = Index
promote_rule(::Type{Poly{A}}, ::Type{Poly{B}}) where {A, B} = Poly{promote_type(A, B)}
promote_rule(::Type{Poly{A}}, ::Type{B}) where {A, B<:Union{ZetaExpr, NN}} = Poly{promote_type(A, B)}

# --- Conversions to Hoffman/Index ---
convert(::Type{Hoffman}, w::HoffmanWord) = Hoffman(Dict{HoffmanWord, Rational{BigInt}}(w => 1))
convert(::Type{Hoffman}, n::NN) = Hoffman(Dict{HoffmanWord, Rational{BigInt}}(HoffmanWord() => convert(Rational{BigInt}, n)))
convert(::Type{Hoffman}, h::Hoffman) = h

convert(::Type{Index}, w::IndexWord) = Index(Dict{IndexWord, Rational{BigInt}}(w => 1))
convert(::Type{Index}, n::NN) = Index(Dict(IndexWord() => convert(Rational{BigInt}, n)))
convert(::Type{Index}, i::Index) = i

# --- Conversions to Poly ---
convert(::Type{Poly{A}}, p::Poly{A}) where A = p
convert(::Type{Poly{A}}, x::Union{ZetaExpr, NN}) where A = Poly{A}(Dict{Int, A}(0 => convert(A, x)))

function convert(::Type{Poly{Hoffman}}, p::Poly{<:Union{NN, Rational{BigInt}}})
    terms = Dict{Int, Hoffman}(d => convert(Hoffman, c) for (d, c) in p.terms)
    return Poly{Hoffman}(terms)
end

function convert(::Type{Poly{Index}}, p::Poly{<:Union{NN, Rational{BigInt}}})
    terms = Dict{Int, Index}(d => convert(Index, c) for (d, c) in p.terms)
    return Poly{Index}(terms)
end

function convert(::Type{Poly{Hoffman}}, p::Poly{Index}; orientation::Symbol = :left)
    terms = Dict{Int, Hoffman}(d => Hoffman(c; orientation=orientation) for (d, c) in p.terms)
    return Poly{Hoffman}(terms)
end

function convert(::Type{Poly{Index}}, p::Poly{Hoffman}; orientation::Symbol = :left)
    terms = Dict{Int, Index}(d => Index(c; orientation=orientation) for (d, c) in p.terms)
    return Poly{Index}(terms)
end

# Generic Poly constructor
function Poly(x::A) where A
    # This is a bit tricky without typelift. We establish a default promotion.
    if A <: Union{Hoffman, HoffmanWord}
        return convert(Poly{Hoffman}, convert(Hoffman, x))
    elseif A <: Union{Index, IndexWord}
        return convert(Poly{Index}, convert(Index, x))
    elseif A <: NN
        return convert(Poly{Rational{BigInt}}, convert(Rational{BigInt}, x))
    else
        return Poly{A}(Dict{Int, A}(0 => x))
    end
end

###################################################################################################
############## Explicit Poly Constructors #########################################################

# --- Poly{Hoffman} ---
Poly{Hoffman}(p::Poly{Hoffman}) = p
Poly{Hoffman}(p::Poly{Index}; orientation::Symbol=:left) = convert(Poly{Hoffman}, p; orientation=orientation)
Poly{Hoffman}(p::Poly{<:Union{NN, Rational{BigInt}}}) = convert(Poly{Hoffman}, p)

Poly{Hoffman}(x::Union{Hoffman, HoffmanWord, NN}) = convert(Poly{Hoffman}, x)
Poly{Hoffman}(x::Union{Index, IndexWord}; orientation::Symbol=:left) = Poly{Hoffman}(Dict(0 => Hoffman(x; orientation=orientation)))

# --- Poly{Index} ---
Poly{Index}(p::Poly{Index}) = p
Poly{Index}(p::Poly{Hoffman}; orientation::Symbol=:left) = convert(Poly{Index}, p; orientation=orientation)
Poly{Index}(p::Poly{<:Union{NN, Rational{BigInt}}}) = convert(Poly{Index}, p)

Poly{Index}(x::Union{Index, IndexWord, NN}) = convert(Poly{Index}, x)
Poly{Index}(x::Union{Hoffman, HoffmanWord}; orientation::Symbol=:left) = Poly{Index}(Dict(0 => Index(x; orientation=orientation)))
