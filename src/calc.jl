#[ calc.jl ]#

# This file implements algorithms for simplifying linear relations of Index objects.
# In particular, it provides elimination procedures that rewrite higher-depth monomial indices
# as linear combinations of lower-depth ones at fixed weight.


#=
export 
=#

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

# function compositions(x::Int)::Vector{Vector{Int}}

# end

"""
    dell_h_1(w::Word) -> Index

```julia
dell_h_1(w) ≡ Index(dell(HoffmanWordtoHoffman(w),1))
```

"""
@inline function dell_h_1(hw::Word)::Index
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
    if get_index_orientation()
        dx = Hoffman([[2,1]])
    else
        dx = Hoffman([[1,2]])
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

#=
-yx&\,y&\,x&\,y&\,x&\,x&\,y&\,x&\,x&\,x \to - &[2, 2, 3, 4]    \\
y&\,-yx&\,x&\,y&\,x&\,x&\,y&\,x&\,x&\,x \to - &[1, 3, 3, 4]    \\
y&\,y&\, yx&\,y&\,x&\,x&\,y&\,x&\,x&\,x \to   &[1, 1, 2, 3, 4] \\
y&\,y&\,x&\,-yx&\,x&\,x&\,y&\,x&\,x&\,x \to - &[1, 2, 4, 4]    \\
y&\,y&\,x&\,y&\, yx&\,x&\,y&\,x&\,x&\,x \to   &[1, 2, 1, 3, 4] \\
y&\,y&\,x&\,y&\,x&\, yx&\,y&\,x&\,x&\,x \to   &[1, 2, 2, 2, 4] \\
y&\,y&\,x&\,y&\,x&\,x&\,-yx&\,x&\,x&\,x \to - &[1, 2, 3, 5]    \\
y&\,y&\,x&\,y&\,x&\,x&\,y&\, yx&\,x&\,x \to   &[1, 2, 3, 1, 4] \\
y&\,y&\,x&\,y&\,x&\,x&\,y&\,x&\, yx&\,x \to   &[1, 2, 3, 2, 3] \\
y&\,y&\,x&\,y&\,x&\,x&\,y&\,x&\,x&\, yx \to   &[1, 2, 3, 3, 2] \\

 yyxyxxyxxx
-y*x*y*x*y*x*x*y*x*x*x \to - &[2, 2, 3, 4]
-y*y*x*x*y*x*x*y*x*x*x \to - &[1, 3, 3, 4]
+y*y*y*x*y*x*x*y*x*x*x \to   &[1, 1, 2, 3, 4]
-y*y*x*y*x*x*x*y*x*x*x \to - &[1, 2, 4, 4]
+y*y*x*y*y*x*x*y*x*x*x \to   &[1, 2, 1, 3, 4]
+y*y*x*y*x*y*x*y*x*x*x \to   &[1, 2, 2, 2, 4]
-y*y*x*y*x*x*y*x*x*x*x \to - &[1, 2, 3, 5]
+y*y*x*y*x*x*y*y*x*x*x \to   &[1, 2, 3, 1, 4]
+y*y*x*y*x*x*y*x*y*x*x \to   &[1, 2, 3, 2, 3]
+y*y*x*y*x*x*y*x*x*y*x \to   &[1, 2, 3, 3, 2]

-y*x*y*x*y*x*x*y*x*x*x-y*y*x*x*y*x*x*y*x*x*x+y*y*y*x*y*x*x*y*x*x*x-y*y*x*y*x*x*x*y*x*x*x+y*y*x*y*y*x*x*y*x*x*x+y*y*x*y*x*y*x*y*x*x*x-y*y*x*y*x*x*y*x*x*x*x+y*y*x*y*x*x*y*y*x*x*x+y*y*x*y*x*x*y*x*y*x*x+y*y*x*y*x*x*y*x*x*y*x

[-y*x*y*x*y*x*x*y*x*x*x,-y*y*x*x*y*x*x*y*x*x*x, y*y*y*x*y*x*x*y*x*x*x,-y*y*x*y*x*x*x*y*x*x*x, y*y*x*y*y*x*x*y*x*x*x, y*y*x*y*x*y*x*y*x*x*x,-y*y*x*y*x*x*y*x*x*x*x, y*y*x*y*x*x*y*y*x*x*x, y*y*x*y*x*x*y*x*y*x*x, y*y*x*y*x*x*y*x*x*y*x]

Index([1,2,3])⨝ Index([5,6,7])

- (Index([1,2,2])⨝ Index([5,6,7]))*up^1
- (Index([1,2,2])⨝ Index([5,6,6]))*up^2
- (Index([1,2,2])⨝ Index([5,6,5]))*up^3
- (Index([1,2,2])⨝ Index([5,6,4]))*up^4
- (Index([1,2,2])⨝ Index([5,6,3]))*up^5
- (Index([1,2,2])⨝ Index([5,6,2]))*up^6
- (Index([1,2,2])⨝ Index([5,6,1]))*up^7
- (Index([1,2,3])⨝ Index([5,6]))*right*up^6



Index([1,2,3])⨝ Index([5,6,6])

- (Index([1,2,2])⨝ Index([5,6,6]))*up^1
- (Index([1,2,2])⨝ Index([5,6,5]))*up^2
- (Index([1,2,2])⨝ Index([5,6,4]))*up^3
- (Index([1,2,2])⨝ Index([5,6,3]))*up^4
- (Index([1,2,2])⨝ Index([5,6,2]))*up^5
- (Index([1,2,2])⨝ Index([5,6,1]))*up^6
- (Index([1,2,3])⨝ Index([5,6]))*right*up^5
=#
