#[ types.jl ]#

# This file contains the type declarations for each

#=
export AbstractOp, OpUp, OpDown, OpLeft, OpRight, OpMinus, OpTau, OpEta, OpPhi, OpDeriv, Operator,
       MPL, ShuffleExpr, HarmonicExpr, ZetaExpr, ExprInt, NN, Word, Hoffman, MonoIndex, Index,
       ShuffleForm, HarmonicForm, MPLCombination, RegHoffman
=#


"""
###################################################################################################
                                        Operators
###################################################################################################
"""


#########################################################
########## Operator System ##############################
# 
# AbstractOp (abstract)
#  ├─ OpUp
#  ├─ OpDown
#  ├─ OpLeft
#  ├─ OpRight
#  ├─ OpMinus
#  ├─ OpTau
#  ├─ OpEta
#  ├─ OpPhi
#  └─ OpDeriv
# 
#########################################################

# 抽象シンボル
abstract type AbstractOp end

mutable struct OpUp    <: AbstractOp
    cnt::Int
end   # ↑
mutable struct OpDown  <: AbstractOp
    cnt::Int
end   # ↓
mutable struct OpLeft  <: AbstractOp
    cnt::Int
end   # ←
mutable struct OpRight <: AbstractOp
    cnt::Int
end   # →
mutable struct OpMinus <: AbstractOp
    cnt::Int
end   # -
mutable struct OpTau   <: AbstractOp
    cnt::Int
end   # τ
mutable struct OpEta   <: AbstractOp
    cnt::Int
end   # η
mutable struct OpPhi   <: AbstractOp
    cnt::Int
end   # φ

# ∂n はパラメータ付き
mutable struct OpDeriv <: AbstractOp
    n::Int
    cnt::Int
end
OpDeriv(n::Int) = OpDeriv(1,n)

function (::Type{T})() where T <: AbstractOp
    if T === OpDeriv
        return T(1,1)
    else
        return T(1)
    end
end

# Operator
mutable struct Operator
    ops::Vector{AbstractOp}

    function Operator()
        new(Vector{AbstractOp}())
    end
end






"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


#########################################################
######## Multiple Polylogarithm System ##################
# 
# MPL (abstract)
#  ├─ ShuffleExpr (abstract)
#  │    └─ ShuffleForm
#  │
#  ├─ HarmonicExpr (abstract)
#  │    ├─ ZetaExpr (abstract)
#  │    │    ├─ Hoffman
#  │    │    ├─ Index
#  │    │    └─ MonoIndex
#  │    │       (Word)
#  │    └─ HarmonicForm
#  │
#  └─ MPLCombination   # 有理数係数上の線形結合
# 
#########################################################


# 大枠：MPL 的なもの（抽象）
abstract type MPL end

# 代数構造で分類（抽象）
abstract type ShuffleExpr  <: MPL end      # 反復積分系（shuffle 構造）
abstract type HarmonicExpr <: MPL end      # z_i=1 に特化（調和和）

# ζ 系（抽象）— 調和和の中でも ζ を中心に
abstract type ZetaExpr     <: HarmonicExpr end

# 指数の型（できたら UInt32/UInt16 へ）
const ExprInt = Int
const NN = Union{Integer,Rational}

# ワード = インデックス列。ハッシュ性と軽さ重視で Tuple に
# WordのTupleの中身はimmutableなものにしておいてください！！！
#const Word = Tuple{Vararg{ExprInt}}  # 例: (2,3) など
# Wordの定義を上記から下記へ変更した このためにWordがTupleのようにふるまうインターフェースをbasefunctions.jlに書いた
struct Word
    t::Tuple{Vararg{Int}}
    function Word(w::Tuple{Vararg{Int}})
        new(w)
    end
end

# 個人によるIndexの向きの補正
const _INDEX_ORIENTATION = Base.RefValue{Bool}(true)
# true  : z_k = y*x^{k-1}
# false : z_k = x^{k-1}*y

set_index_orientation!(val::Bool) = (_INDEX_ORIENTATION[] = val)
get_index_orientation() = _INDEX_ORIENTATION[]

# Hoffman 代数の元：ワードの有限線形結合（係数は有理数）
# xy^3x^2 -> [1,2,2,2,1,1]
mutable struct Hoffman <: ZetaExpr
    terms::Dict{Word,Rational{BigInt}}
    function Hoffman()
        new(Dict{Word,Rational{BigInt}}())
    end
end

# 「ζ系の生のインデックス」を薄いラッパで持つ
mutable struct MonoIndex <: ZetaExpr
    word::Word
    coeff::Rational{BigInt}
    function MonoIndex(word::Word, coeff::NN)
        if iszero(coeff)
            throw(DomainError(coeff,"coefficient cant not be 0."))
        end
        new(word,coeff)
    end
end

# yxyx^3y^2x^2 -> [2,4,1,3]
# もし _INDEX_ORIENTATION がfalseなら
# x^2y^2x^3yxy -> [3,1,4,2]
mutable struct Index <: ZetaExpr
    terms::Dict{Word, Rational{BigInt}}
    function Index()
        new(Dict{Word, Rational{BigInt}}())
    end
end

# 一般の shuffle 形式
struct ShuffleForm <: ShuffleExpr

end
# 一般の harmonic 形式
struct HarmonicForm <: HarmonicExpr

end

# “MPL1つが複数の表現の有理係数線形結合”を許す総和器
struct MPLCombination <: MPL
    terms::Dict{MPL,Rational{BigInt}}     # 異なる具象（Hoffman/EulerSum/ShuffleForm）も混在OK
end

# Hoffman型の値を係数とする多項式
mutable struct RegHoffman
    terms::Dict{Int,Hoffman}  # T のべき -> Hoffman 元の係数
    function RegHoffman()
        new(Dict{Int,Hoffman}())
    end
end
