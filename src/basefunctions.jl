#[ basefunctions.jl ]#

# This file defines the basic functions
import Base: iszero, isone, zero, one, copy, ==, lastindex
import Base: getindex, length, iterate, firstindex, eltype, isempty, copy, collect, Tuple, vcat, hcat, reverse, ==, hash, isless, isequal, convert, in, iterate, haskey, keys, values, pairs, keytype, valtype

#=
export is_monomial, is_hoffmanword, is_hoffman, is_index, is_indexword,
       is_shuffleform, is_harmonicform, is_mplcombination, is_shuffleexpr, is_harmonicexpr, is_zetaexpr,
       is_admissible
=#

"""
###################################################################################################
                                        Arithmatics
###################################################################################################
"""
#=
A
C(n+1,k+1) = (n+1)/(k+1)*C(n,k)
B
C(n+1,k) = (n+1)/(n-k+1)*C(n,k)
=#
#= (遅い)
function multinomial(v::Vector{Int})
    sorted_v = sort(v)
    binomm = Rational(BigInt(1))
    n = sorted_v[1]
    k = sorted_v[1]
    
    multinomm = Rational(BigInt(1))

    for i in 1:lastindex(v)-1
        for _ in 1:(sorted_v[i+1]-sorted_v[i])
            n += 1
            k += 1
            binomm *= n//k
        end
        for _ in 1:sorted_v[i]
            n += 1
            binomm *= n//(n-k)
        end
        multinomm *= binomm
    end

    return multinomm
end
=#
@inline function multinomial(v::Vector{Int})::BigInt
    num = factorial(BigInt(sum(v)))
    den = BigInt(1)
    for k in v
        den *= factorial(BigInt(k))
    end
    return div(num,den)
end
"""
v,n -> n/(v[1]!v[2]! ...)
"""
@inline function multinomial(v::Vector{Int},n::BigInt)::BigInt
    den = BigInt(1)
    for k in v
        den *= factorial(BigInt(k))
    end
    return div(n,den)
end


"""
###################################################################################################
                                        Operators
###################################################################################################
"""

###################################################################################################
############## Property Functions #################################################################


isone(op::Operator)::Bool = isempty(op.ops)
isone(op::Ts where Ts<:AbstractOp)::Bool  = (op.cnt == 0)

one(::Type{Operator})::Operator = Operator()
one(op::Type{<:AbstractOp})::AbstractOp = (op)(1)

==(a::Operator, b::Operator)::Bool = a.ops == b.ops
function ==(a::AbstractOp, b::AbstractOp)::Bool
    ta = typeof(a)
    tb = typeof(b)
    if ta === tb
        if ta === OpDeriv
            if a.cnt == b.cnt && a.n == b.n
                return true
            end
        else
            if a.cnt == b.cnt
                return true
            end
        end
        return false
    elseif (ta === OpUp && tb === OpDown) || (ta === OpDown && tb == OpUp)
        if a.cnt == -b.cnt
            return true
        end
    end
    return false
end
==(a::Operator, b::AbstractOp)::Bool = lastindex(a.ops) == 1 && a.ops[1] == b

lastindex(op::Operator)::Int = lastindex(op.ops)

function ==(a::Hoffman, b::NN)
    if iszero(a)
        return iszero(b)
    end
    if is_monomial(a)
        (k, v) = first(a.terms)
        return isone(k) && v == b
    end
    return false
end
==(a::NN, b::Hoffman) = b == a

function ==(a::Index, b::NN)
    if iszero(a)
        return iszero(b)
    end
    if is_monomial(a)
        (k, v) = first(a.terms)
        return isone(k) && v == b
    end
    return false
end
==(a::NN, b::Index) = b == a

###################################################################################################
############## Base Functions #####################################################################

# 文字->Operatorに変換するテーブル
const _String_to_Operator_Table = Dict{String, DataType}(
    "↑" => OpUp,
    "↓" => OpDown,
    "←" => OpLeft,
    "→" => OpRight,
    "-" => OpMinus,
    "τ" => OpTau,
    "†" => OpTau, # constで使おうとすると unknown unicode character '†'
    "𖼷" => OpTau, # ミャオ文字なので使える(泣)
    "η" => OpEta,
    "⋁" => OpEta,
    "φ" => OpPhi,
    "ϕ" => OpPhi,
    "∂" => OpDeriv,
)
# Operator->文字に変換するテーブル
const _Operator_to_String_Table = Dict{DataType, String}(
    OpUp => "↑",
    OpDown => "↓",
    OpLeft => "←",
    OpRight => "→",
    OpMinus => "-",
    OpTau => "τ",
    OpEta => "η",
    OpPhi => "φ",
    OpDeriv => "∂",
)

const _RE_DERIV = r"^∂(\d+)$"

@inline function str_to_op(sym::AbstractString,n::Int = 1)::AbstractOp
    m = match(_RE_DERIV, sym)
    if m !== nothing
        i = parse(Int,m.captures[1])
        return OpDeriv(n,i)
    end
    return _String_to_Operator_Table[sym](n)
end
@inline function str_to_op(t::Tuple{AbstractString,Integer})::AbstractOp
    m = match(_RE_DERIV, t[1])
    if m !== nothing
        i = parse(Int,m.captures[1])
        return OpDeriv(t[2],i)
    end
    return _String_to_Operator_Table[t[1]](t[2])
end

@inline function clean(ops_in::Vector{<:AbstractOp})::Vector{AbstractOp}
    isempty(ops_in) && return AbstractOp[]

    # 1. Normalize all ops (e.g., OpDown -> OpUp)
    ops = [copy(op) for op in ops_in]

    # 2. Merge adjacent ops
    cleaned_ops = AbstractOp[]
    sizehint!(cleaned_ops, length(ops))

    for current_op in ops
        if isempty(cleaned_ops)
            push!(cleaned_ops, current_op)
            continue
        end

        last_op = cleaned_ops[end]
        T_last = typeof(last_op)
        T_curr = typeof(current_op)

        merged = if T_last !== T_curr
            # Different types, cannot merge
            nothing
        elseif T_curr === OpTau
            # xor for OpTau
            OpTau((last_op.cnt + current_op.cnt) % 2)
        elseif T_curr === OpDeriv && last_op.n == current_op.n
            # Merge derivatives of same order
            OpDeriv(last_op.cnt + current_op.cnt, last_op.n)
        elseif T_curr !== OpDeriv
            # Merge other ops by adding counts
            T_curr(last_op.cnt + current_op.cnt)
        else
            # Different order derivatives, cannot merge
            nothing
        end

        if merged !== nothing
            # Replace the last element with the merged result
            cleaned_ops[end] = merged
        else
            # Append the new, unmergable op
            push!(cleaned_ops, current_op)
        end

        # If the last op's count became zero, remove it
        if !isempty(cleaned_ops) && cleaned_ops[end].cnt == 0
            pop!(cleaned_ops)
        end
    end
    return cleaned_ops
end

@inline function copy(op::Operator)::Operator
    r = Operator()
    if isone(op)
        return r
    end

    # ops is a vector of immutable structs, so a shallow copy is sufficient.
    return Operator(copy(op.ops))
end
# これは正規化も兼ねている copy
@inline function copy(op::AbstractOp)::AbstractOp
    top = typeof(op)
    if top === OpDeriv
        return OpDeriv(op.cnt,op.n)
    elseif top === OpDown
        return OpUp(-op.cnt)
    elseif top === OpTau
        return OpTau(op.cnt % 2)
    else
        return (top)(op.cnt)
    end
end




"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## Word compatible interface ##########################################################


# ===== Tuple 互換 =====
getindex(w::AbstractWord, i::Int) = w.t[i]
getindex(w::AbstractWord, r::UnitRange) = (typeof(w))(w.t[r])
length(w::AbstractWord) = length(w.t)
iterate(w::AbstractWord, s...) = iterate(w.t, s...)
firstindex(w::AbstractWord) = firstindex(w.t)
lastindex(w::AbstractWord) = lastindex(w.t)
eltype(::Type{AbstractWord}) = Int
eltype(::Type{HoffmanWord})  = Int
eltype(::Type{IndexWord})    = Int
isempty(w::AbstractWord) = isempty(w.t)

# ===== 各種操作互換 =====
@inline copy(w::W) where {W<:AbstractWord} = W(w.t)                 # コピー
collect(w::AbstractWord) = collect(Int,w.t)                         # Vector化
Tuple(w::AbstractWord) = w.t                                        # タプル化
vcat(a::W, b::W) where {W<:AbstractWord} = W((a.t..., b.t...))      # 連結
reverse(w::W) where {W<:AbstractWord} = W(reverse(w.t))             # 反転
==(a::W, b::W) where {W<:AbstractWord} = a.t == b.t
hash(w::AbstractWord, h::UInt) = hash(w.t, h)
isless(a::W, b::W) where {W<:AbstractWord} = isless(a.t,b.t)
isequal(a::W, b::W) where {W<:AbstractWord} = isequal(a.t, b.t)
function convert(::Type{W}, t::Tuple{Vararg{Int}})::W where {W<:AbstractWord}
    return W(t)
end
# ===== スプラット互換 (a... が動くように) =====
in(item, w::AbstractWord) = in(item, w.t)
iterate(w::AbstractWord) = iterate(w.t)  # ←これが超重要！


###################################################################################################
############## Hoffman compatible interface ##########################################################

length(a::Hoffman) = length(a.terms)
length(a::Index)   = length(a.terms)
length(a::Poly)    = length(a.terms)

isempty(a::Hoffman) = isempty(a.terms)
isempty(a::Index)   = isempty(a.terms)
isempty(a::Poly)    = isempty(a.terms)

keys(a::Hoffman) = keys(a.terms)
keys(a::Index)   = keys(a.terms)
keys(a::Poly)    = keys(a.terms)

keytype(::Type{Hoffman})         = HoffmanWord
keytype(::Type{Index})           = IndexWord
keytype(::Type{Poly{A}}) where A = Int

values(a::Hoffman) = values(a.terms)
values(a::Index)   = values(a.terms)
values(a::Poly)    = values(a.terms)

valtype(::Type{Hoffman})         = Rational{BigInt}
valtype(::Type{Index})           = Rational{BigInt}
valtype(::Type{Poly{A}}) where A = A


pairs(a::Hoffman) = pairs(a.terms)
pairs(a::Index)   = pairs(a.terms)
pairs(a::Poly)    = pairs(a.terms)

iterate(a::Hoffman) = iterate(a.terms)
iterate(a::Index)   = iterate(a.terms)
iterate(a::Poly)    = iterate(a.terms)

iterate(a::Hoffman,s...) = iterate(a.terms,s...)
iterate(a::Index,s...)   = iterate(a.terms,s...)
iterate(a::Poly,s...)    = iterate(a.terms,s...)

eltype(::Type{Hoffman})         = Tuple{HoffmanWord,Rational{BigInt}}
eltype(::Type{Index})           = Tuple{IndexWord  ,Rational{BigInt}}
eltype(::Type{Poly{A}}) where A = Tuple{Int,A}

getindex(h::Hoffman, w::HoffmanWord) = get(h.terms, w, zero(Rational{BigInt}))
getindex(i::Index, w::IndexWord)     = get(i.terms, w, zero(Rational{BigInt}))
getindex(r::Poly{A}, d::Int) where A = get(r.terms, d, zero(A))

haskey(h::Hoffman, w::HoffmanWord) = haskey(h.terms, w)
haskey(i::Index, w::IndexWord)     = haskey(i.terms, w)
haskey(r::Poly, d::Int)            = haskey(r.terms, d)

==(a::Hoffman, b::Hoffman)         = a.terms == b.terms
==(a::Index, b::Index)             = a.terms == b.terms
==(a::Poly{A}, b::Poly{A}) where A = a.terms == b.terms

isequal(a::Hoffman, b::Hoffman)         = isequal(a.terms, b.terms)
isequal(a::Index, b::Index)             = isequal(a.terms, b.terms)
isequal(a::Poly{A}, b::Poly{A}) where A = isequal(a.terms, b.terms)

hash(h::Hoffman, u::UInt) = hash(h.terms, u)
hash(i::Index, u::UInt)   = hash(i.terms, u)
hash(r::Poly, u::UInt)    = hash(r.terms, u)

convert(::Type{Hoffman}, d::Dict{HoffmanWord, Rational{BigInt}})::Hoffman = Hoffman(d)
convert(::Type{Index}, d::Dict{IndexWord, Rational{BigInt}})::Index       = Index(d)
function convert(::Type{Poly{A}}, d::Dict{Int, A})::Poly{A} where A
    return Poly{A}(d)
end

@inline copy(a::Hoffman)::Hoffman = Hoffman(Dict(a.terms))
@inline copy(a::Index)::Index = Index(Dict(a.terms))
@inline function copy(p::Poly{A})::Poly{A} where A
    r = Poly{A}()
    for (deg,h) in p
        r.terms[deg] = copy(h)
    end
    return r
end


###################################################################################################
############## Property Functions #################################################################


iszero(w::AbstractWord)::Bool     = false
iszero(i::Index)::Bool            = isempty(i.terms)
iszero(w::Hoffman)::Bool          = isempty(w.terms)
#iszero(hrm::HarmonicForm)::Bool   =  # TODO
#iszero(shf::ShuffleForm)::Bool    =  # TODO
#iszero(mpl::MPLCombination)::Bool =  # TODO
iszero(r::Poly)::Bool             = isempty(r.terms)

isone(w::AbstractWord)::Bool     = isempty(w.t)
isone(w::Hoffman)::Bool          = length(w) == 1 && haskey(w,HoffmanWord()) && isone(w[HoffmanWord()])
isone(i::Index)::Bool            = length(i) == 1 && haskey(i,IndexWord()) && isone(i[IndexWord()])
#isone(hrm::HarmonicForm)::Bool   =  # TODO
#isone(shf::ShuffleForm)::Bool    =  # TODO
#isone(mpl::MPLCombination)::Bool =  # TODO
isone(r::Poly)::Bool             = length(r) == 1 && haskey(r,0) && isone(r[0])

function is_admissible(idx::Index; orientation::Symbol = :left)::Bool
    if orientation === :left # z_k = y*x^{k-1}
        for w in keys(idx)
            if !isempty(w) && w[end] <= 1
                return false
            end
        end
        return true
    elseif orientation === :right # z_k = x^{k-1}*y
        for w in keys(idx)
            if !isempty(w) && w[1] <= 1
                return false
            end
        end
        return true
    end
    throw(ArgumentError("orientation must be :left or :right"))
end
function is_admissible(h::Hoffman; orientation::Symbol = :left)::Bool
    if orientation === :left # Admissible words must end with x0
        for w in keys(h)
            # Admissible words are of the form y*...*x, which means they start with 1 and end with 0
            if !isempty(w) && (w[1] != 1 || w[end] != 0)
                return false
            end
        end
        return true
    elseif orientation === :right # Admissible words must start with x0
        for w in keys(h)
            # Admissible words are of the form x*...*y, which means they start with 0 and end with 1
            if !isempty(w) && (w[1] != 0 || w[end] != 1)
                return false
            end
        end
        return true
    end
    throw(ArgumentError("orientation must be :left or :right"))
end

zero(::Type{Index})::Index               = Index()
zero(::Type{Hoffman})::Hoffman           = Hoffman()
zero(::Type{HarmonicForm})::HarmonicForm = HarmonicForm()
zero(::Type{ShuffleForm})::ShuffleForm   = ShuffleForm()
function zero(::Type{Poly{A}})::Poly{A} where A
    return Poly{A}()
end
one(::Type{HoffmanWord})::HoffmanWord = HoffmanWord()
one(::Type{IndexWord})::IndexWord     = IndexWord()
function one(::Type{Index})::Index
    idx = Index()
    idx.terms[IndexWord()] = Rational(BigInt(1))
    return idx
end
function one(::Type{Hoffman})::Hoffman
    w = Hoffman()
    w.terms[HoffmanWord()] = Rational(BigInt(1))
    return w
end
function one(::Type{Poly{A}})::Poly{A} where A
    r = Poly{A}()
    r.terms[0] = one(A)
    return r
end
# Instance methods
zero(::Index) = zero(Index)
zero(::Hoffman) = zero(Hoffman)
zero(::HarmonicForm) = zero(HarmonicForm)
zero(::ShuffleForm) = zero(ShuffleForm)
zero(::Poly{A}) where A = zero(Poly{A})

one(::HoffmanWord) = one(HoffmanWord)
one(::IndexWord) = one(IndexWord)
one(::Index) = one(Index)
one(::Hoffman) = one(Hoffman)
one(::Poly{A}) where A = one(Poly{A})
# TODO: one for hrm shf mpl

is_monomial(i::Index)::Bool        = length(i) == 1
is_monomial(w::Hoffman)::Bool      = length(w) == 1
is_monomial(w::AbstractWord)::Bool = true
#is_monomial(hrm::HarmonicForm)::Bool        = # TODO
#is_monomial(shf::ShuffleForm)::Bool         = # TODO
#is_monomial(mpl::MPLCombination)::Bool      = # TODO
is_monomial(r::Poly)::Bool         = length(r) == 1

is_hoffmanword(a)::Bool    = typeof(a) === HoffmanWord
is_hoffman(a)::Bool        = typeof(a) === Hoffman
is_indexword(a)::Bool      = typeof(a) === IndexWord
is_index(a)::Bool          = typeof(a) === Index
is_shuffleform(a)::Bool    = typeof(a) === ShuffleForm
is_harmonicform(a)::Bool   = typeof(a) === HarmonicForm
is_mplcombination(a)::Bool = typeof(a) === MPLCombination
is_shuffleexpr(a)::Bool    = typeof(a) <:  ShuffleExpr
is_harmonicexpr(a)::Bool   = typeof(a) <:  HarmonicExpr
is_zetaexpr(a)::Bool       = typeof(a) <:  ZetaExpr


###################################################################################################
############## Base Functions #####################################################################

# index compressor
""" example (1,1,0,0,1,0) -> [1,3,2] """
@inline function idxprs(w::HoffmanWord)::Vector{Int}
    n2 = count(==(1), w)
    v = Vector{Int}(undef,n2)
    idx = 0
    for x in w
        if x == 1
            idx += 1
            v[idx] = 1
        else
            v[idx] += 1
        end
    end
    return v
end
""" example (0,1,0,0,1,1) -> [2,3,1] """
@inline function idxprs_r(w::HoffmanWord)::Vector{Int}
    n2 = count(==(1), w)
    v = Vector{Int}(undef,n2)
    idx = 0
    cnt = 1
    for x in w
        if x == 0
            cnt += 1
        else
            idx += 1
            v[idx] = cnt
            cnt = 1
        end
    end
    return v
end
# index expander
""" example [1,3,2] -> (1,1,0,0,1,0) """
@inline function idxdprs(v::Vector{Int})::HoffmanWord
    # 総長さを一度に計算してから確保
    total_len = sum(v)
    w = zeros(Int,total_len)

    pos = 1
    for t in v
        w[pos] = 1
        pos += t
    end
    return HoffmanWord(Tuple(w))
end
""" example [2,3,1] -> (0,1,0,0,1,1) """
@inline function idxdprs_r(v::Vector{Int})::HoffmanWord
    # 総長さを一度に計算してから確保
    total_len = sum(v)
    w = zeros(Int,total_len)

    pos = 0
    for t in v
        pos += t
        w[pos] = 1
    end
    return HoffmanWord(Tuple(w))
end

# Helper to add terms to a dictionary
function _add!(d::Dict{W, C}, w::W, c::C) where {W, C}
    iszero(c) && return
    new_val = get(d, w, zero(C)) + c
    if iszero(new_val)
        delete!(d, w)
    else
        d[w] = new_val
    end
end

@inline function _add_terms(a_terms::Dict{K,V}, b_terms::Dict{K,V}, op) where {K,V}
    new_terms = copy(a_terms)
    for (k, v) in b_terms
        new_val = op(get(new_terms, k, zero(V)), v)
        if iszero(new_val)
            delete!(new_terms, k)
        else
            new_terms[k] = new_val
        end
    end
    return new_terms
end