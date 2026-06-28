#[ mathcore.jl ]#

# This file defines the fundamental functions for numerical computation.



function multinomial(v::Vector{Int})::BigInt
    num = factorial(BigInt(sum(v)))
    den = BigInt(1)
    for k in v
        den *= factorial(BigInt(k))
    end
    return div(num,den)
end
"""
v -> sum(v)!/(v[1]!v[2]! ...)

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
    partition_count(n::Int)
五角数定理を用いて分割数 p(n) を計算する。sizehint! 用。
"""
function partition_count(n::Int)::Int
    n < 0 && return 0
    n == 0 && return 1
    p = zeros(Int, n+1)
    p[1] = 1
    for i in 2:n+1
        k1 = 0
        k2 = 0
        k_diff = 1
        i_sign = true
        while true
            k1 += k_diff
            k2 += k_diff + 1
            k_diff += 3
            if i < k2 + 1
                break
            end
            if i_sign
                p[i] += p[i - k1]
                p[i] += p[i - k2]
                i_sign = false
            else
                p[i] -= p[i - k1]
                p[i] -= p[i - k2]
                i_sign = true
            end
        end
        if i > k1
            if i_sign
                p[i] += p[i - k1]
            else
                p[i] -= p[i - k1]
            end
        end
    end
    return p[n+1]
end

"""
    partitions(n::Int)
n の分割をすべて返す。Zoghbi and Stojmenovic アルゴリズムの改良版。
"""
function partitions(n::Int)::Vector{Vector{Int}}
    n < 0 && return Vector{Vector{Int}}()
    n == 0 && return [Int[]]
    
    result = Vector{Vector{Int}}()

    # 事前にサイズを予約
    sizehint!(result, partition_count(n))
    
    # 内部バッファ
    ps = ones(Int, n)
    ps[1] = n
    m = 1 # 現在の分割の長さ
    h = 1 # 1より大きい要素の右端インデックス
    
    # 最初の分割 [n] を追加
    push!(result, ps[1:m])
    
    while ps[1] != 1
        if ps[h] == 2
            ps[h] = 1
            h -= 1
            m += 1
        else
            r = ps[h] - 1
            t = m - h + 1
            ps[h] = r
            
            # 残りの値を r で埋め尽くす
            while t >= r
                h += 1
                ps[h] = r
                t -= r
            end
            
            if t == 0
                m = h
            else
                m = h + 1
                if t > 1
                    h += 1
                    ps[h] = t
                end
            end
        end
        # ps[1:m] は Julia ではコピーを作成するので、これで安全かつ高速
        push!(result, ps[1:m])
    end
    
    return result
end


function compositions(n::Int, k::Int; init::Union{Nothing,Vector{Int}}=nothing)::Vector{Vector{Int}}

    if k < 1 || n < k
        return Vector{Vector{Int}}()
    end

    if k == 1
        v = [n]
        if init !== nothing
            v[1] += init[1]
        end
        return [v]
    end

    len = binomial(n-1,k-1) 
    results = Vector{Vector{Int}}(undef,len)

    curr = ones(Int,k)
    curr[k]=n-k+1

    ###################
    # initなし（最速）
    ###################

    if init === nothing

        results[1]=copy(curr)

        for i in 2:len
            j=k-1
            while j>0 && curr[j+1]==1
                j-=1
            end
            curr[j]+=1
            excess=curr[j+1]-1
            curr[j+1]=1
            curr[k]=excess
            results[i]=copy(curr)

        end

        return results

    end

    ###################
    # initあり
    ###################

    @inbounds begin
        v = copy(curr)
        for i=1:k
            v[i]+=init[i]
        end
        results[1]=v

        for r in 2:len

            j=k-1
            while j>0 && curr[j+1]==1
                j-=1
            end

            curr[j]+=1
            excess=curr[j+1]-1
            curr[j+1]=1
            curr[k]=excess
            v = copy(curr)

            for i=1:k
                v[i]+=init[i]
            end
            results[r]=v
        end
    end
    return results
end

function compositions(n::Int)::Vector{Vector{Int}}
    n > 0 || throw(ArgumentError("n must be positive"))
    
    # 2^(n-1) 個の要素を事前に確保
    num_compositions = 1 << (n - 1)
    result = Vector{Vector{Int}}(undef, num_compositions)
    
    for i in 0:(num_compositions - 1)
        comp = Int[]
        # 最悪ケース（1+1+...+1）を考慮して容量を確保
        sizehint!(comp, n)

        current_sum = 1
        # i の bit 単位で仕切りを判定
        for j in 0:(n - 2)
            if (i >> j) & 1 == 1
                push!(comp, current_sum)
                current_sum = 1
            else
                current_sum += 1
            end
        end
        push!(comp, current_sum)
        result[i + 1] = comp
    end
    
    return result
end



"""
Zagier予想に基づく重さ k の MZV 空間の次元 d_k を求める。
d_k = d_{k-2} + d_{k-3}, d_0=1, d_1=0, d_2=1
"""
function zagier_dim(k::Int)
    k < 0 && return 0
    k == 0 && return BigInt(1)
    k == 1 && return BigInt(0)
    k == 2 && return BigInt(1)

    # 空間計算量 O(1): 直近3つの値のみを保持
    # d_p3: d_{i-3}, d_p2: d_{i-2}, d_p1: d_{i-1}
    d_p3 = BigInt(1) # d_0
    d_p2 = BigInt(0) # d_1
    d_p1 = BigInt(1) # d_2

    for _ in 3:k
        # d_i = d_{i-2} + d_{i-3}
        d_curr = d_p2 + d_p3
        
        # 参照の更新（シフト）
        d_p3 = d_p2
        d_p2 = d_p1
        d_p1 = d_curr
    end
    
    return d_p1
end

