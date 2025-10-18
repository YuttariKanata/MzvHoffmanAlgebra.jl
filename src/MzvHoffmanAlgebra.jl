module MzvHoffmanAlgebra


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




       # types.jl
export AbstractOp, OpUp, OpDown, OpLeft, OpRight, OpMinus, OpTau, OpEta, OpPhi, OpDeriv, Operator,
       MPL, ShuffleExpr, HarmonicExpr, ZetaExpr, ExprInt, NN, Word, Hoffman, MonoIndex, Index,
       ShuffleForm, HarmonicForm, MPLCombination, RegHoffman,
       # basefunctions.jl
       is_monomial, is_monoindex, is_hoffman, is_index,
       is_shuffleform, is_harmonicform, is_mplcombination, is_shuffleexpr, is_harmonicexpr,is_zetaexpr,
       # converting.jl
       HoffmanWordtoMonoIndex, IndexWordtoMonoIndex,
       word, IndexWordtoHoffmanWord, HoffmanWordtoIndexWord,
       HoffmanWordtoIndex, IndexWordtoIndex,
       HoffmanWordtoHoffman, IndexWordtoHoffman,
       x, y, T,
       # arithmetic.jl
       shift_degree, add!,
       # hoffman.jl
       monomial_sh, monomial_st, monomial_sh_double, monomial_st_double,
       shuffle_product, stuffle_product, shuffle_product_double, stuffle_product_double,
       shuffle_pow, stuffle_pow, shpw,
       Hoffman_hom, Hoffman_antihom, monomial_sw_w, starword_to_word,
       monomial_dual, dual,
       # accessors.jl
       upper_represent, sortedprint,
       # operator.jl
       left_act, right_act, ⬆️, ➡️, ⬇️, ⬅️, ➖, up, right, down, left, minus, τ, 𖼷, η, ⋁, φ, ∂

include("types.jl")
include("basefunctions.jl")
include("converting.jl")
include("arithmetic.jl")
include("hoffman.jl")
include("accessors.jl")
include("operator.jl")

end
