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
       MPL, ShuffleExpr, HarmonicExpr, ZetaExpr, NN, AbstractWord, HoffmanWord, Hoffman, IndexWord, Index,
       ShuffleForm, HarmonicForm, MPLCombination, Poly, T

       # basefunctions.jl
export multinomial, str_to_op, clean,
       is_monomial, is_hoffmanword, is_hoffman, is_index, is_indexword,
       is_shuffleform, is_harmonicform, is_mplcombination, is_shuffleexpr, is_harmonicexpr, is_zetaexpr,
       is_admissible,
       idxprs, idxprs_r, idxdprs, idxdprs_r

       # converting.jl
export index, x, y

       # arithmetic.jl
export shift_degree

       # products.jl
export shuffle_product, stuffle_product, star_stuffle_product,
       ш, ⨝, ∗, ⋆,
       shuffle_pow, stuffle_pow, star_stuffle_pow,
       st_index1_pow, sh_y_pow

       # regularization.jl
export stuffle_regularization, shuffle_regularization, stuffle_regularization_polynomial, shuffle_regularization_polynomial, rho, rho_t

       # duals.jl
export dual, Hoffman_dual, Landen_dual,
       Hoffman_hom, Hoffman_antihom, starword_to_word

       # derivations.jl
export dell, Hoffman_derivation

       # accessors.jl
export sortedprint, upper_represent

       # operators.jl
export left_act, right_act, ⬆️, ➡️, ⬇️, ⬅️, ➖, up, right, down, left, minus, τ, 𖼷, η, ⋁, φ, ∂,
       WordtoOperator

       # interporation.jl
export shift, sharp_shift
;


include("types.jl")
include("basefunctions.jl")
include("converting.jl")
include("arithmetic.jl")
include("monomials.jl")
include("products.jl")
include("regularization.jl")
include("duals.jl")
include("derivations.jl")
include("accessors.jl")
include("operators.jl")
include("interporation.jl")

end