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
       ShuffleForm, HarmonicForm, MPLCombination, Poly, T,
       set_index_orientation!, get_index_orientation,
       # basefunctions.jl
       is_monomial, is_hoffman, is_index, is_monoindex,
       is_shuffleform, is_harmonicform, is_mplcombination, is_shuffleexpr, is_harmonicexpr,is_zetaexpr,
       # converting.jl
       HoffmanWordtoMonoIndex, IndexWordtoMonoIndex,
       IndexWordtoHoffmanWord, HoffmanWordtoIndexWord,
       HoffmanWordtoIndex, IndexWordtoIndex,
       HoffmanWordtoHoffman, IndexWordtoHoffman,
       x, y,
       # arithmetic.jl
       shift_degree, add!,
       # hoffman.jl
       shuffle_product, stuffle_product, star_stuffle_product, 
       shuffle_product_double, stuffle_product_double, star_stuffle_product_double,
       shuffle_pow, stuffle_pow, star_stuffle_pow, 
       shpw, stpw, starstpw,
       Hoffman_hom, Hoffman_antihom, starword_to_word,
       dual, Hoffman_dual, Landen_dual,
       # accessors.jl
       upper_represent, sortedprint,
       # operator.jl
       left_act, right_act, ⬆️, ➡️, ⬇️, ⬅️, ➖, up, right, down, left, minus, τ, 𖼷, η, ⋁, φ, ∂,
       WordtoOperator

include("types.jl")
include("basefunctions.jl")
include("converting.jl")
include("arithmetic.jl")
include("hoffman.jl")
include("accessors.jl")
include("operator.jl")

end
