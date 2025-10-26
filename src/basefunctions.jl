#[ basefunctions.jl ]#

# This file defines the basic functions
import Base: iszero, isone, zero, one, copy, ==, lastindex
import Base: getindex, length, iterate, firstindex, eltype, isempty, copy, collect, Tuple, vcat, hcat, reverse, ==, hash, one, iterate, show, isless, isequal, in, convert

#=
export is_monomial, is_monoindex, is_hoffman, is_index,
       is_shuffleform, is_harmonicform, is_mplcombination, is_shuffleexpr, is_harmonicexpr,is_zetaexpr
=#

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
            if a.n == b.n && a.cnt == b.cnt
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

###################################################################################################
############## Base Functions #####################################################################


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

@inline function str_to_op(sym::AbstractString,n::Int = 1)::AbstractOp
    m = match(r"^∂(\d+)$", sym)
    if m !== nothing
        i = parse(Int,m.captures[1])
        return OpDeriv(i,n)
    end
    return _String_to_Operator_Table[sym](n)
end
@inline function str_to_op(t::Tuple{AbstractString,Integer})::AbstractOp
    m = match(r"^∂(\d+)$", t[1])
    if m !== nothing
        i = parse(Int,m.captures[1])
        return OpDeriv(i,t[2])
    end
    return _String_to_Operator_Table[t[1]](t[2])
end


@inline function append_clean!(op::Operator, vn::Vector{<:AbstractOp})::Operator

    vnl = lastindex(vn)
    v = Vector{AbstractOp}(undef,vnl)
    for i in 1:vnl
        v[i] = copy(vn[i])
    end

    # 常に op.ops内のUp,DownはすべてUpに直す
    for vi in v

        bef = (isempty(op.ops) ? nothing : op.ops[end])

        tvi = typeof(vi)
        top = typeof(bef)
        if tvi !== top
            push!(op.ops, vi)
        elseif tvi === OpTau
            op.ops[end].cnt = xor(vi.cnt, bef.cnt) & 1
        elseif tvi === OpDeriv
            if vi.n == bef.n
                op.ops[end].cnt += vi.cnt
            else
                push!(op.ops,vi)
            end
        else
            op.ops[end].cnt += vi.cnt
        end

        if !isempty(op.ops) && op.ops[end].cnt == 0
            pop!(op.ops)
        end

    end

    return op
end
@inline function append_clean!(a::Operator, bn::Operator)::Operator
    
    aftend = lastindex(bn.ops)
    b = Vector{AbstractOp}(undef,aftend)
    for i in 1:aftend
        x = bn.ops[i]
        if x isa OpDeriv
            b[i] = OpDeriv(x.n,x.cnt)
        else
            b[i] = typeof(x)(x.cnt)
        end
    end
    if isempty(b)
        return a
    end
    if isempty(a.ops)
        r = Operator()
        r.ops = b
        return r
    end

    # 中央境界
    befidx = lastindex(a.ops)
    aftidx = 1

    bef = a.ops[befidx]
    aft = b[aftidx]

    aftend += 1

    while true
        # 全てUpに正規化されている
        tb = typeof(bef)
        ta = typeof(aft)
        if tb !== ta
            break  # 中央で融合不可能なら終了
        elseif tb === OpTau
            bef.cnt = xor(bef.cnt, aft.cnt) & 1
        elseif tb === OpDeriv
            if bef.n == aft.n
                bef.cnt += aft.cnt
            end
        else
            bef.cnt += aft.cnt
        end

        # もしcnt==0なら消滅 → 再帰的処理へ
        if bef.cnt == 0
            befidx -= 1
            aftidx += 1
            if aftidx == aftend
                break
            end
            if befidx == 0
                break
            end
            bef = a.ops[befidx]
            aft = b[aftidx]
            continue
        else
            # 結合成功 → b側先頭を削除し終了
            a.ops[befidx] = bef
            aftidx += 1
            break
        end
    end

    resize!(a.ops,befidx)
    # 残りを連結
    append!(a.ops,@view b[aftidx:aftend-1])
    return a
end

@inline function copy(op::Operator)::Operator
    r = Operator()
    if isone(op)
        return Operator()
    end

    opl = lastindex(op.ops)
    b = Vector{AbstractOp}(undef,opl)
    for i in 1:opl
        x = op.ops[i]
        if x isa OpDeriv
            b[i] = OpDeriv(x.n,x.cnt)
        else
            b[i] = typeof(x)(x.cnt)
        end
    end
    append!(r.ops,b)
    return r
end
# これは正規化も兼ねている copy
@inline function copy(op::AbstractOp)::AbstractOp
    top = typeof(op)
    if top === OpDeriv
        return OpDeriv(op.n,op.cnt)
    elseif top === OpDown
        return OpUp(-op.cnt)
    elseif top === OpTau
        return OpTau(op.cnt & 1)
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
getindex(w::Word, i::Int) = w.t[i]
getindex(w::Word, r::UnitRange) = Word(w.t[r])
length(w::Word) = length(w.t)
iterate(w::Word, s...) = iterate(w.t, s...)
firstindex(w::Word) = firstindex(w.t)
lastindex(w::Word) = lastindex(w.t)
eltype(::Type{Word}) = ExprInt
isempty(w::Word) = isempty(w.t)

# ===== 各種操作互換 =====
copy(w::Word) = Word(w.t)                               # コピー
collect(w::Word) = collect(w.t)                         # Vector化
Tuple(w::Word) = w.t                                   # タプル化
vcat(a::Word, b::Word) = Word((a.t..., b.t...))         # 連結
hcat(a::Word, b::Word) = Word((a.t..., b.t...))
reverse(w::Word) = Word(reverse(w.t))                   # 反転
==(a::Word, b::Word) = a.t == b.t
hash(w::Word, h::UInt) = hash(w.t, h)
isless(a::Word, b::Word) = isless(a.t,b.t)
isequal(a::Word, b::Word) = isequal(a.t, b.t)
convert(::Type{Word}, t::Tuple)::Word = Word(t)

# ===== スプラット互換 (a... が動くように) =====
iterate(w::Word) = iterate(w.t)  # ←これが超重要！
in(item, w::Word) = in(item, w.t)



###################################################################################################
############## Property Functions #################################################################


iszero(w::Word)::Bool             = false
iszero(x::Index)::Bool            = isempty(x.terms)
iszero(m::MonoIndex)::Bool        = false
iszero(w::Hoffman)::Bool          = isempty(w.terms)
#iszero(hrm::HarmonicForm)::Bool   =  # TODO
#iszero(shf::ShuffleForm)::Bool    =  # TODO
#iszero(mpl::MPLCombination)::Bool =  # TODO
iszero(r::RegHoffman)::Bool       = isempty(r.terms)

isone(w::Word)::Bool             = w == Word()
isone(x::Index)::Bool            = length(x.terms) == 1 && haskey(x.terms,Word()) && isone(x.terms[Word()])
isone(m::MonoIndex)::Bool        = isone(m.word) && isone(m.coeff)
isone(w::Hoffman)::Bool          = length(w.terms) == 1 && haskey(w.terms,Word()) && isone(w.terms[Word()])
#isone(hrm::HarmonicForm)::Bool   =  # TODO
#isone(shf::ShuffleForm)::Bool    =  # TODO
#isone(mpl::MPLCombination)::Bool =  # TODO
isone(r::RegHoffman)::Bool       = length(r.terms) == 1 && haskey(r.terms,0) && isone(r.terms[0])

==(a::Index,b::Index)::Bool           = a.terms == b.terms
==(a::MonoIndex,b::MonoIndex)::Bool   = a.word == b.word && a.coeff == b.coeff
==(a::Hoffman,b::Hoffman)::Bool       = a.terms == b.terms
==(a::RegHoffman,b::RegHoffman)::Bool = a.terms == b.terms

zero(::Type{T}) where T <: Union{Index,Hoffman,HarmonicForm,ShuffleForm,MPLCombination} = T()
one(::Type{Word})::Word = Word()
one(::Type{MonoIndex})::MonoIndex = MonoIndex(Word(),1)
function one(::Type{Index})::Index
    idx = Index()
    idx.terms[Word()] = 1
    return  idx
end
function one(::Type{Hoffman})::Hoffman
    w = Hoffman()
    w.terms[Word()] = 1
    return  w
end
function one(::Type{RegHoffman})::RegHoffman
    r = RegHoffman()
    r.terms[0] = one(Hoffman)
    return r
end
# TODO: one for hrm shf mpl

is_monomial(x::Index)::Bool                 = isone(length(x.terms))
is_monomial(w::Hoffman)::Bool               = isone(length(w.terms))
is_monomial(w::Union{Word,MonoIndex})::Bool = true
#is_monomial(hrm::HarmonicForm)::Bool      = # TODO
#is_monomial(shf::ShuffleForm)::Bool       = # TODO
#is_monomial(mpl::MPLCombination)::Bool    = # TODO
is_monoindex(r::RegHoffman)::Bool           = isone(length(r.terms))

is_hoffman(x::MPL)::Bool        = typeof(x) == Hoffman
is_index(x::MPL)::Bool          = typeof(x) == Index
is_monoindex(x::MPL)::Bool      = typeof(x) == MonoIndex
is_shuffleform(x::MPL)::Bool    = typeof(x) == ShuffleForm
is_harmonicform(x::MPL)::Bool   = typeof(x) == HarmonicForm
is_mplcombination(x::MPL)::Bool = typeof(x) == MPLCombination
is_shuffleexpr(x::MPL)::Bool    = typeof(x) <: ShuffleExpr
is_harmonicexpr(x::MPL)::Bool   = typeof(x) <: HarmonicExpr
is_zetaexpr(x::MPL)::Bool       = typeof(x) <: ZetaExpr


###################################################################################################
############## Base Functions #####################################################################

# index compressor
""" example (2,2,1,1,2,1) -> [1,3,2] """
@inline function idxprs(w::Word)::Vector{Int}
    n2 = count(==(2), w)
    v = Vector{Int}(undef,n2)
    idx = 0
    for x in w
        if x == 2
            idx += 1
            v[idx] = 1
        else
            v[idx] += 1
        end
    end
    return v
end
# index expander
""" example [1,3,2] -> (2,2,1,1,2,1) """
@inline function idxdprs(v::Vector{Int})::Word
    # 総長さを一度に計算してから確保
    total_len = sum(v)
    w = Vector{ExprInt}(undef, total_len)

    pos = 1
    for t in v
        w[pos] = 2
        pos += 1
        if t > 0
            fill!(view(w, pos:pos+t-2), 1)
            pos += t-1
        end
    end
    return Word(w)
end

@inline function copy(a::Hoffman)::Hoffman
    r = Hoffman()
    for (w,c) in a.terms
        r.terms[w] = c
    end
    return r
end
@inline function copy(a::Index)::Index
    r = Index()
    for (w,c) in a.terms
        r.terms[w] = c
    end
    return r
end
@inline function copy(p::RegHoffman)::RegHoffman
    r = RegHoffman()
    for (deg,h) in p.terms
        r.terms[deg] = copy(h)
    end
    return r
end