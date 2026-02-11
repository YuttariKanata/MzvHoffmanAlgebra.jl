#[ types.jl ]#

# This file contains the type declarations for each

#=
export AbstractOp, OpUp, OpDown, OpLeft, OpRight, OpMinus, OpTau, OpEta, OpPhi, OpDeriv, Operator,
       MPL, ShuffleExpr, HarmonicExpr, ZetaExpr, NN, AbstractWord, HoffmanWord, Hoffman, IndexWord, Index,
       ShuffleForm, HarmonicForm, MPLCombination, Poly, T,
       set_index_orientation!, get_index_orientation
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

# cntは繰り返しの回数
struct OpUp    <: AbstractOp
    cnt::Int
end   # ↑
struct OpDown  <: AbstractOp
    cnt::Int
end   # ↓
struct OpLeft  <: AbstractOp
    cnt::Int
end   # ←
struct OpRight <: AbstractOp
    cnt::Int
end   # →
struct OpMinus <: AbstractOp
    cnt::Int
end   # -
struct OpTau   <: AbstractOp
    cnt::Int
end   # τ
struct OpEta   <: AbstractOp
    cnt::Int
end   # η
struct OpPhi   <: AbstractOp
    cnt::Int
end   # φ

# ∂n はパラメータ付き
struct OpDeriv <: AbstractOp
    cnt::Int
    n::Int
end
OpDeriv(n::Int) = OpDeriv(n,1)

(::Type{A})() where {A<:AbstractOp} = A(1)

# Operator
struct Operator
    ops::Vector{AbstractOp}
    function Operator(v::Vector{<:AbstractOp})
        new(clean(v))
    end
end

Operator() = Operator(AbstractOp[])




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
#  │    └─ HarmonicForm
#  │
#  └─ MPLCombination   # 有理数係数上の線形結合
#
# ZetaExpr (abstract)
#   ├─ Hoffman
#   ├─ Index
#   ├─ MonoIndex
#   └─ Word 
#########################################################


# 大枠：MPL 的なもの（抽象）
abstract type MPL end

# 代数構造で分類（抽象）
abstract type ShuffleExpr  <: MPL end      # 反復積分系（shuffle 構造）
abstract type HarmonicExpr <: MPL end      # z_i=1 に特化（調和和）

# ζ 系（抽象）— 調和和の中でも ζ を中心に
abstract type ZetaExpr end

# 有理数と整数
const NN = Union{Integer,Rational}

# ワード = インデックス列。ハッシュ性と軽さ重視で Tuple に
# WordのTupleの中身はimmutableなものにしておいてください！！！
#const Word = Tuple{Vararg{Int}}  # 例: (2,3) など
# Wordの定義を上記から下記へ変更した このためにWordがTupleのようにふるまうインターフェースをbasefunctions.jlに書いた
abstract type AbstractWord <: ZetaExpr end

struct HoffmanWord <: AbstractWord
    t::Tuple{Vararg{Int}}
end
struct IndexWord <: AbstractWord
    t::Tuple{Vararg{Int}}
end

HoffmanWord() = HoffmanWord(())
IndexWord()   = IndexWord(())

# The orientation of indices (e.g. :left for z_k=y*x^{k-1} or :right for z_k=x^{k-1}*y)
# is now passed as an argument to conversion functions, instead of using a global state.

# REPLで表示する項の数(目安)
const _OMIT_COUNTS = 100

# Hoffman 代数の元：ワードの有限線形結合（係数は有理数）
# xy^3x^2 -> [1,2,2,2,1,1]
struct Hoffman <: ZetaExpr
    terms::Dict{HoffmanWord,Rational{BigInt}}
    function Hoffman(d::Dict{HoffmanWord,Rational{BigInt}})
        new(filter(p->!iszero(p.second),d))
    end
end
Hoffman() = Hoffman(Dict{HoffmanWord,Rational{BigInt}}())

# yxyx^3y^2x^2 -> [2,4,1,3]
# もし _INDEX_ORIENTATION がfalseなら
# x^2y^2x^3yxy -> [3,1,4,2]
struct Index <: ZetaExpr
    terms::Dict{IndexWord, Rational{BigInt}}
    function Index(d::Dict{IndexWord,Rational{BigInt}})
        new(filter(p->!iszero(p.second),d))
    end
end
Index() = Index(Dict{IndexWord,Rational{BigInt}}())

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

# 多項式一つで済ます
struct Poly{A}
    terms::Dict{Int,A}
    function Poly{A}(d::Dict{Int,A}) where A
        new{A}(filter(p->!iszero(p.second),d))
    end
end
Poly{A}() where A = Poly{A}( Dict{Int,A}() )

const T::Poly{Rational{BigInt}} = Poly{Rational{BigInt}}( Dict{Int,Rational{BigInt}}( 1 => Rational{BigInt}(1) ) )