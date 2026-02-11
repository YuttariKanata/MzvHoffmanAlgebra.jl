#[ calc.jl ]#

# This file implements algorithms for simplifying linear relations of Index objects.
# In particular, it provides elimination procedures that rewrite higher-depth monomial indices
# as linear combinations of lower-depth ones at fixed weight.


#=
export 
=#

"""
       ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
       ┃ This file is under development. ┃
       ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
"""


function compositions(x::Int)::Vector{Vector{Int}}
    x > 0 || throw(ArgumentError("x must be positive"))

    result = Vector{Vector{Int64}}()

    function dfs(remaining::Int, current::Vector{Int})
        if remaining == 0
            push!(result, copy(current))
            return
        end
        for k in 1:remaining
            push!(current, k)
            dfs(remaining - k, current)
            pop!(current)
        end
    end

    dfs(x, Int[])
    return result
end

"""
    dell_h_1(w::Word) -> Index

```julia
    dell_h_1(w) ≡ Index(dell(Hoffman(w),1))
```

"""
@inline function dell_h_1(hw::HoffmanWord; orientation=:left)::Index
    #=
    s = Hoffman()
    lw = lastindex(w)
    for i in 1:lw
        add!(s,w[1:i-1]* image[w[i]] * w[i+1:lw])
    end
    return s
    =#
    s = Hoffman()
    lw = lastindex(hw)
    if orientation == :left
        dx = Hoffman([[2,1]])
    elseif orientation == :right
        dx = Hoffman([[1,2]])
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end
    for i in 1:lw
        if hw[i] == 1
            #@show hw[1:i-1], dx, hw[i+1:lw]
            add!(s,  hw[1:i-1] * dx * hw[i+1:lw] )
        else
            #@show hw[1:i-1], -dx, hw[i+1:lw]
            add!(s,-(hw[1:i-1] * dx * hw[i+1:lw]))
        end
    end
    return Index(s)
end

"""
    admissible_indices(n::Int, k::Int)::Vector{Vector{Int}}

weight = n, depth = k の許容インデックス（先頭 ≥ 2）をすべて返す
"""
function admissible_indices(n::Int, k::Int)::Vector{Vector{Int}}
    if k <= 0 || n < k + 1
        return Vector{Vector{Int}}()
    end

    # 初期値 [2,1,1,...,1]
    a = ones(Int64,k)
    a[1] = 2

    r = n - (k + 1)   # 余剰

    # 結果サイズは binomial(n-2, k-1)
    res = Vector{Vector{Int}}()
    sizehint!(res, binomial(BigInt(n - 2), k - 1))

    while true
        # 出力（コピーはここだけ）

        # 右から余剰を動かす
        i = k
        while i >= 1 && r == 0
            r += a[i] - 1
            a[i] = 1
            i -= 1
        end

        if i < 1
            break
        end

        a[i] += 1
        r -= 1

        @show a

        push!(res, copy(a))

    end

    return res
end

function bm(u)
    res =BitVector(undef, sizeof(u)*8) 
    res.chunks[1] = u%UInt64
    res
end

function solve_liner_dependencies(n::Int)
    for k in n-1:-1:1
        # インデックスの深さがkのインデックスが現れる関係式を見つけていく
        
        # 導分関係式 dell(w,1) ( =-reg_sh(y∗w) )
        # dell(w,1)が返すインデックスのweightは wt(w)+1 depthは dep(w) or dep(w)+1
        # depthがkの関係式が欲しいなら入力はdepthがk-1のものにする
        # weightはn-1のもの

        

    end
end