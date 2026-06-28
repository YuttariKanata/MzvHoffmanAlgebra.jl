module MzvHoffmanAlgebra

using REPL
using Dates


# モジュールのメインコンテンツ...

function __init__()  # moduleが読み込まれた後に最初に実行される特別な関数
    # REPLがインタラクティブな場合のみ表示（スクリプト実行時は邪魔にならないように）
    if isinteractive() && Base.get(stdout, :color, false)
        print_banner()
    end
end

function print_banner()
    # 81 (Light Cyan) -> 76 (Spring Green)
    start_color = 81
    end_color = 76
    
    # ASCII Art の作成（Hoffmanの 'H' や MZV を意識したデザイン）
    # ここでは数学的な「重み」や「深さ」をイメージした抽象的なデザインにしている
    logo = raw"""
     __  _________    __     _    _        __  __                       
    |  \/ /___\ \/   / /    | |  | |      / _|/ _|                      
    | \  / |  /\ \  / /_____| |__| | ___ | |_| |_ _ __ ___   __ _ _ __ 
    | |\/| | / /\ \/ /______|  __  |/ _ \|  _|  _| '_ ` _ \ / _` | '_ \ 
    | |  | |/ /__\_ /       | |  | | (_) | | | | | | | | | | (_| | | | |
    |_|  |_/______//        |_|  |_|\___/|_| |_| |_| |_| |_|\__,_|_| |_|
    
    """
    
    # 1行ずつ装飾して出力
    lines = split(logo, '\n')
    println()
    for line in lines
        isempty(line) && continue
        len = length(line)
        for (i, char) in enumerate(line)
            # 線形補間で色を決定 (i=1 で 81, i=len で 76)
            # 数値が減る方向なので max(76, 81 - ...) で制御
            t = (i - 1) / (len > 1 ? len - 1 : 1)
            c = round(Int, start_color + t * (end_color - start_color))
            
            printstyled(char, color=c, bold=true)
        end
        println()
    end
    println()
    printstyled("  ", ">>> Multiple Zeta Values & Hoffman Algebra Package\n\n", color=:white)

    printstyled("  GitHub: https://github.com/YuttariKanata/MzvHoffmanAlgebra.jl\n\n", color=248)

    # バージョン情報や更新日の表示（任意）
    printstyled("  $(rpad(string(pkgversion(MzvHoffmanAlgebra)),10)) | $(today()) | Happy Computing with MZVs!\n", color=:light_black)
    println("-"^70)
end



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


include("mathcore.jl")
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
include("calc.jl")


end