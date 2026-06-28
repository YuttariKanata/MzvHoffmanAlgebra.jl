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
        dx = Hoffman([[1,0]])
    elseif orientation == :right
        dx = Hoffman([[0,1]])
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end
    for i in 1:lw
        if hw[i] == 1
            #@show hw[1:i-1], dx, hw[i+1:lw]
            s += hw[1:i-1] * dx * hw[i+1:lw]
        else
            #@show hw[1:i-1], -dx, hw[i+1:lw]
            s -= hw[1:i-1] * dx * hw[i+1:lw]
        end
    end
    return Index(s)
end

"""
    admissible_indices(n::Int)::Vector{Vector{Int}}

weight = n の許容インデックス（最後 ≥ 2）をすべて返す
"""
function admissible_indices(n::Int; orientation::Symbol=:left)::Vector{Vector{Int}}
    # n=1 の場合は空集合（許容インデックスは存在しない）
    n > 1 || return Vector{Vector{Int}}()
    
    # 重さ n-1 の全組成の個数は 2^(n-2)
    num_indices = 1 << (n - 2)
    result = Vector{Vector{Int}}(undef, num_indices)
    if orientation == :left
        for i in 0:(num_indices - 1)
            idx = Int[]
            sizehint!(idx, n)
            
            current_sum = 1  # 先頭は 1 から開始
            
            # 重さ n-1 の組成を作るための (n-2) 回のビット走査
            for j in 0:(n - 3)
                if (i >> j) & 1 == 1
                    push!(idx, current_sum)
                    current_sum = 1
                else
                    current_sum += 1
                end
            end
            
            # 最後の成分に +1 して、全体の重さを n にし、かつ末尾 >= 2 を保証
            push!(idx, current_sum + 1)
            result[i + 1] = idx
        end
    elseif orientation == :right
        for i in 0:(num_indices - 1)
            # 内部的な push! 回数を減らすため sizehint!
            # 最大要素数は (n-1) なので n-1 を確保しておけば十分
            idx = Int[]
            sizehint!(idx, n-1)
            
            # 組成生成ロジック (重さ n-1)
            # current_sum の初期値を 1 ではなく 2 にすることで、
            # 先頭要素に +1 した状態を直接作る
            current_sum = 2
            
            # ビット走査範囲は (n-1) 個の 1 の間の隙間 (n-2) 箇所
            for j in 0:(n - 3)
                if (i >> j) & 1 == 1
                    push!(idx, current_sum)
                    current_sum = 1
                else
                    current_sum += 1
                end
            end
            push!(idx, current_sum)
            result[i + 1] = idx
        end
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end

    return result
end

# 非アロケーション・イテレータ
function admissible_indices(n::Int, k::Int)::Vector{Vector{Int}}
    # 許容条件: 重さ n, 深さ k, 末尾 >= 2
    # これは重さ m = n-1, 深さ k, 末尾 >= 1 の組成に 1 を足すのと同値
    if k < 1 || n < k + 1
        return Vector{Vector{Int}}()
    end

    m = n - 1
    result = Vector{Vector{Int}}()
    
    # 辞書式順序の最初: [m-k+1, 1, 1, ..., 1]
    # (先頭に重さを寄せ、右側を最小の 1 で埋める)
    current = ones(Int, k)
    current[1] = m - k + 1
    
    # 常に末尾以外にある「1より大きい要素」の右端インデックスを追跡してもいいが、
    # シンプルな走査でも十分速い。ここでは完全な非再帰ループで記述。
    while true
        # 1. 現在の組成をコピーして保存 (末尾に +1 して重さ n に調整)
        res = copy(current)
        @inbounds res[k] += 1
        push!(result, res)
        
        # 2. 次の組成への遷移 (辞書式順序)
        # 右から見て「1より大きい」かつ「末尾ではない」要素を探す
        idx = k - 1
        @inbounds while idx > 0 && current[idx] == 1
            idx -= 1
        end
        
        # 全て 1 なら走査完了
        idx == 0 && break
        
        # 3. 遷移処理: idx の値を 1 減らし、その右隣 (idx+1) に残りの重さを全て集約する
        # そのさらに右側 (idx+2 .. k) は 1 に戻す
        @inbounds begin
            current[idx] -= 1
            # idx より右側に割り当てるべき合計重さを計算
            # 遷移前の段階で current[idx+1...k] は 1...1, (余り) だったので、
            # 新しい current[idx+1] は (現在の余り + 1) になる
            rem_weight = current[k] + 1
            current[idx+1] = rem_weight
            
            # idx+1 より右側がある場合は、そこを 1 にリセット
            if idx + 1 < k
                current[k] = 1
            end
        end
    end
    
    return result
end

function bm(u)
    res =BitVector(undef, sizeof(u)*8) 
    res.chunks[1] = u%UInt64
    res
end


f1(x) = shuffle_product(index(1),x) - stuffle_product(index(1),x)
f2(x) = shuffle_product(index(2),x) - stuffle_product(index(2),x)
f3(x) = shuffle_product(index(3),x) - stuffle_product(index(3),x)
f5(x) = shuffle_product(index(4),x) - stuffle_product(index(4),x)
f4(x) = shuffle_product(index(1,2),x) - stuffle_product(index(1,2),x)

function sb(A::AbstractMatrix, pad_width::Int)
    rows, cols = size(A)
    println("size:",size(A))
    print(" "^pad_width)
    for i in 1:size(A,2)
        print("$(lpad(i,pad_width)):")
    end
    println()
    println("_"^((pad_width+1)*(size(A,2)+1)+2))
    for i in 1:rows
        # 行の開始
        
        print("$(lpad(i,pad_width))[")
        
        for j in 1:cols
            val = A[i, j]
            
            # 1. Rational型なら可能ならIntに変換
            if val isa Rational
                if denominator(val) == 1
                    val = numerator(val)
                end
            end
            
            # 2. 0は表示しない、それ以外は文字列化
            str_val = (val == 0) ? "" : string(val)

            str_val *= ":"
            
            # 3. lpadで整列して出力（最後の列以外はスペースを少し空ける）
            print(lpad(str_val, pad_width))
            if j < cols
                print(" ") 
            end
        end
        
        # 行の終了
        println(" ]$i")
    end
end
function sb2(A::AbstractMatrix, b::AbstractVector, pad_width::Int)
    rows, cols = size(A)
    println("size:",size(A))
    for i in 1:rows
        # 行の開始
        print("$(lpad(i,pad_width))[")
        
        for j in 1:cols
            val = A[i, j]
            
            # 1. Rational型なら可能ならIntに変換
            if val isa Rational
                if denominator(val) == 1
                    val = numerator(val)
                end
            end
            
            # 2. 0は表示しない、それ以外は文字列化
            str_val = (val == 0) ? "" : string(val)
            
            # 3. lpadで整列して出力（最後の列以外はスペースを少し空ける）
            print(lpad(str_val, pad_width))
            if j < cols
                print(" ") 
            end
        end
        
        # 行の終了
        print(" ] ")
        println(b[i])
    end
end

function f(n,x::Union{Index,IndexWord})
    # if n==1
    #     index(1) ⨝ x - index(1) ∗ x
    #     # shuffle_product(index(1),x) - stuffle_product(index(1),x)
    # elseif n==2        
    #     index(2) ⨝ x - index(2) ∗ x
    #     # shuffle_product(index(2),x) - stuffle_product(index(2),x)
    # elseif n==3
    #     index(3) ⨝ x - index(3) ∗ x
    #     # shuffle_product(index(3),x) - stuffle_product(index(3),x)
    # elseif n==4
    #     index(4) ⨝ x - index(4) ∗ x
    #     # shuffle_product(index(4),x) - stuffle_product(index(4),x)
    # elseif n==5
    #     index(5) ⨝ x - index(5) ∗ x
    # end
    n ⨝ x - n ∗ x
end

function max_depth(i::Index)
    m = -1
    for (w,c) in i.terms
        m=max(m,length(w))
    end
    return m
end

"""
Hoffman's relation (k_r >= 2 convention):
LHS (depth r) - RHS (depth r+1) = 0
Input: k = (k1, k2, ..., kr) where kr >= 2
"""
function hoffman_relation(k::Vector{Int})
    r = length(k)
    rel = Index() 

    # --- LHS: \sum_{i=1}^r \zeta(k1, ..., ki+1, ..., kr) ---
    # 全ての項で最後尾は kr または kr+1 なので、kr >= 2 なら許容的
    for i in 1:r
        ki_plus_1 = copy(k)
        ki_plus_1[i] += 1
        #add_term!(rel, IndexWord(ki_plus_1), 1)
        rel += IndexWord(ki_plus_1)
    end

    # --- RHS: \sum_{i=1}^r \sum_{j=1}^{ki-1} \zeta(k1, ..., j, ki-j+1, ..., kr) ---
    for i in 1:r
        val_i = k[i]
        # k_i = 1 のときは内側の和は空（j=1:0）
        for j in 1:(val_i - 1)
            # 新しい深さ r+1 のベクトルを pre-allocate
            new_word_vec = Vector{Int}(undef, r + 1)
            
            # 低レベルな詰め込み
            @inbounds begin
                # i-1 番目までコピー
                for l in 1:(i-1)
                    new_word_vec[l] = k[l]
                end
                
                # 分割部分: ki -> (j, ki-j+1)
                new_word_vec[i] = j
                new_word_vec[i+1] = val_i - j + 1
                
                # i+1 番目からコピー
                for l in (i+1):r
                    new_word_vec[l+1] = k[l]
                end
            end
            
            # i=r の場合、最後尾は k_r - j + 1
            # j <= k_r - 1 なので、k_r - j + 1 >= 2 が保証され、admissible となる
            # add_term!(rel, IndexWord(new_word_vec), -1)
            rel -= IndexWord(new_word_vec)
        end
    end

    return rel
end

#=

Index(1,1,1,7) + 1/6*Index(3,7) + 1/24*Index(10) - 1/2*Index(1,7,2) + 1/4*Index(2,2,6) + 1/24*Index(4,6) - 1/8*Index(6,4) + 1/2*Index(2,1,7) + 1/2*Index(1,1,8) + 1/6*Index(1,9) + 1/2*Index(1,2,1,6) - 1/2*Index(1,1,6,2) + 1/4*Index(6,2,2) - 1/3*Index(7,3) - 1/4*Index(2,6,2) + 1/4*Index(2,8) + 1/2*Index(1,1,2,6) + 1/6*Index(3,1,6) + 1/6*Index(1,3,6) + Index(1,1,1,1,6) - 1/3*Index(1,6,3) + 1/2*Index(2,1,1,6) + 1/2*Index(1,2,7) - 1/4*Index(8,2)

=#

function solve_mzv2(n::Int)
    rels = Vector{Index}()
    admis = admissible_indices(n)
    sort!(admis, by = x -> (length(x), x), rev = true)
    d=Dict{IndexWord,Int}()
    for i in eachindex(admis)
        d[index(admis[i])] = i
    end
    
    function dell_reduction(x::Index)
        r = copy(x)
        for (w,c) in x.terms
            if d[w] <= BigInt(2)^(n-3)
                r -= w*c
                r += dual(w)*c
            end
        end
        return r
    end
    
    cnt = 0
    for ad_is in reverse(admissible_indices(n-1))
        

        rel = dell_reduction(f(index(1),index(ad_is)))
        if rel != 0
            cnt += 1
            if cnt % 2 == 0
                continue
            end
            push!(rels,rel)
        end
        #push!(rels,(f(index(1),index(ad_is))))
    end

    #sort!(rels, by = x->(-max_depth(x),minimum([d[idx] for idx in keys(x.terms)])), rev = false)


    for i in 2:3
        for w1 in admissible_indices(i), w2 in admissible_indices(n-i)
            rel = dell_reduction(f(index(w1),index(w2)))
            if rel != 0
                push!(rels,rel)
            end
            #push!(rels,f(index(w1),index(w2)))
        end
    end
    # for i in 1:n
    #     for ad_is in reverse(admissible_indices(n-i))
    #         push!(rels,f(index(i),index(ad_is)))
    #     end
    # end
    for ad_is in admissible_indices(n)
        rel = Index(ad_is)-dual(Index(ad_is))
        if !isnothing(rel) && rel != 0
            push!(rels,Index(ad_is)-dual(Index(ad_is)))
        end
    end

    sort!(rels, by = x->(-max_depth(x),minimum([d[idx] for idx in keys(x.terms)])), rev = false)



    return rels
end


function solve_mzv_basis(n::Int)
    # 1. 全ての許容インデックスを取得
    admis = admissible_indices(n)
    
    # 2. 列の順序付け: 深さが「大きい」ものを左（消去対象）、
    # 「小さい」ものを右（基底候補）に配置する。
    # 深さが同じなら辞書式順序で並べる。
    sort!(admis, by = x -> (length(x), x), rev = true)
    
    admisdict = Dict{IndexWord, Int}(IndexWord(admis[i]) => i for i in eachindex(admis))
    ncols = length(admis)
    
    # 3. 関係式の収集 (君の solve_mzv2 を利用)
    # 注: hoffman_relation は weight n-1 -> n のものを想定
    # solve_mzv2(n::Int)::Vector{Index} Indexの和が0になるものを集めている
    raw_rels = solve_mzv2(n) 
    
    # 4. 行列の構築
    # A * x = 0 の形式。行数は関係式の数、列数は admis の数
    nrows = length(raw_rels)
    A = zeros(Rational{BigInt}, nrows, ncols)
    
    for (r_idx, rel) in enumerate(raw_rels)
        for (word, coeff) in rel
            c_idx = get(admisdict, word, nothing)
            if !isnothing(c_idx)
                A[r_idx, c_idx] += coeff
            end
        end
    end


    # 5. ガウスの消去法 (Row Echelon Form)
    # ピボット選択を行い、上三角に近い形へ
    pivots = Vector{Int}() # 各行のピボット列インデックス
    curr_row = 1
    for j in 1:ncols
        #sb(A,9)
        # j列目で非ゼロを持つ行を curr_row 以降から探す
        #pivot_row = A[curr_row, j] != 0 ? 1 : nothing
        pivot_row = findfirst(i -> A[i, j] != 0, curr_row:nrows)
        #@show pivot_row A[curr_row, j]
        if isnothing(pivot_row)
            continue # この列は自由変数（基底の一部）
        end
        pivot_row += (curr_row - 1)
        
        # 行入れ替え
        A[curr_row, :], A[pivot_row, :] = A[pivot_row, :], A[curr_row, :]
        push!(pivots, j)
        
        # 正規化 (A[curr_row, j] を 1 に)
        pivot_val = A[curr_row, j]
        A[curr_row, :] .//= pivot_val
        
        # 他の行の掃き出し (Reduced Row Echelon Form を目指す)
        for i in 1:nrows
            if i != curr_row && A[i, j] != 0
                factor = A[i, j]
                A[i, :] -= factor * A[curr_row, :]
            end
        end
        prev_pivot = pivot_val
        curr_row += 1
        curr_row > nrows && break
    end


    # 6. 結果の整理
    # ピボットが立った列 = 他のインデックスで書けるもの
    # ピボットが立たなかった列 = 基底 (Basis)
    basis_indices = setdiff(1:ncols, pivots)
    
    solutions = Dict{IndexWord, Index}()
    for (idx, p_col) in enumerate(pivots)
        # s_{p_col} = - \sum_{non_pivot} A[idx, j] * s_j
        expr = Index()
        for j in (p_col + 1):ncols
            if A[idx, j] != 0
                #add_term!(expr, IndexWord(admis[j]), -A[idx, j])
                expr -= copy( IndexWord(admis[j]) * A[idx, j] )
            end
        end
        solutions[IndexWord(admis[p_col])] = expr
    end

    return solutions, [admis[b] for b in basis_indices]
end

"""
大野の関係式を生成する。
k: 重さ n-m の許容インデックス (Vector{Int})
m: 増加させる重さ
"""
function ohno_relation(k::Vector{Int}, m::Int)
    r = length(k)
    k_dual = collect(dual(index(k))) # Vector{Int} を返すと想定
    r_dual = length(k_dual)
    
    rel = Index()
    
    # LHS: Σ ζ(k + e)
    for e in multiset_partitions(m, r)
        word = k .+ e
        # kr >= 2 かつ er >= 0 ならば word[end] >= 2 は維持される
        #add_term!(rel, IndexWord(word), 1)
        rel += IndexWord(word)
    end
    
    # RHS: Σ ζ(k_dual + e')
    for e_prime in multiset_partitions(m, r_dual)
        word_dual = k_dual .+ e_prime
        #add_term!(rel, IndexWord(word_dual), -1)
        rel -= (IndexWord(word_dual))
    end
    
    return rel
end

function multiset_partitions(n::Int, k::Int)
    results = Vector{Vector{Int}}()
    current = zeros(Int, k)
    
    function backtrack(remain, idx)
        if idx == k
            current[idx] = remain
            push!(results, copy(current))
            return
        end
        for val in 0:remain
            current[idx] = val
            backtrack(remain - val, idx + 1)
        end
    end
    
    backtrack(n, 1)
    return results
end

#=

    print(mzv(7) - mzv(2,5) - mzv(5,2) - mzv(1,6) - mzv(3,4) - mzv(4,3))
    print(mzv(2,5) - mzv(1,1,5) + mzv(1,6) - mzv(1,2,4) - mzv(1,3,3) - mzv(1,4,2))
    print(-mzv(2,3,2) - mzv(2,1,4) + mzv(2,5) + mzv(3,4) - mzv(1,2,4) - mzv(2,2,3))
    print(mzv(2,1,4) + mzv(1,1,5) - mzv(1,1,1,4) - mzv(1,1,2,3) + mzv(1,2,4) - mzv(1,1,3,2))
    print(-mzv(3,2,2) - mzv(3,1,3) + mzv(4,3) + mzv(3,4) - mzv(1,3,3) - mzv(2,2,3))
    print(-mzv(1,2,1,3) - mzv(1,1,2,3) - mzv(1,2,2,2) + mzv(1,2,4) + mzv(1,3,3) + mzv(2,2,3))
    print(-mzv(1,2,1,3) - mzv(2,1,2,2) + mzv(2,1,4) + mzv(3,1,3) - mzv(2,1,1,3) + mzv(2,2,3))
    print(mzv(1,2,1,3) - mzv(1,1,1,1,3) + mzv(1,1,1,4) + mzv(2,1,1,3) + mzv(1,1,2,3) - mzv(1,1,1,2,2))
    print(-mzv(2,3,2) + mzv(5,2) - mzv(3,2,2) + mzv(4,3) - mzv(4,1,2) - mzv(1,4,2))
    print(mzv(2,3,2) - mzv(1,3,1,2) - mzv(1,2,2,2) + mzv(1,3,3) - mzv(1,1,3,2) + mzv(1,4,2))
    print(mzv(2,3,2) - mzv(2,1,2,2) + mzv(3,2,2) - mzv(1,2,2,2) - mzv(2,2,1,2) + mzv(2,2,3))
    print(mzv(2,1,2,2) - mzv(1,1,2,1,2) + mzv(1,1,2,3) + mzv(1,2,2,2) - mzv(1,1,1,2,2) + mzv(1,1,3,2))
    print(mzv(3,2,2) - mzv(3,1,1,2) + mzv(3,1,3) - mzv(1,3,1,2) - mzv(2,2,1,2) + mzv(4,1,2))
    print(mzv(1,2,1,3) - mzv(1,1,2,1,2) + mzv(1,3,1,2) + mzv(2,2,1,2) + mzv(1,2,2,2) - mzv(1,2,1,1,2))
    print(mzv(2,1,2,2) - mzv(2,1,1,1,2) + mzv(3,1,1,2) + mzv(2,1,1,3) + mzv(2,2,1,2) - mzv(1,2,1,1,2))
    print(mzv(2,1,1,1,2) + mzv(1,1,2,1,2) + mzv(1,1,1,1,3) + mzv(1,2,1,1,2) + mzv(1,1,1,2,2) - mzv(1,1,1,1,1,2))
    print(3*mzv(2,4) - mzv(6) + 2*mzv(3,3) + 8*mzv(1,5))
    print(-mzv(3,3) + 3*mzv(1,2,3) + 9*mzv(1,1,4) - mzv(1,5))
    print(-mzv(2,4) - mzv(4,2) + 4*mzv(1,2,3) + 4*mzv(1,3,2) + 4*mzv(2,1,3))
    print(8*mzv(1,1,1,3) + mzv(1,2,1,2) - mzv(3,1,2) - mzv(1,3,2) + 2*mzv(1,1,2,2) - mzv(1,1,4))
    print(6*mzv(2,4) - mzv(6) + 12*mzv(1,5))
    print(-mzv(4,2) + 4*mzv(1,2,3) + mzv(1,3,2) + mzv(2,2,2) + 9*mzv(1,1,4) - mzv(1,5) + 2*mzv(2,1,3))
    print(-mzv(4,2) + 4*mzv(1,2,3) + mzv(1,3,2) + mzv(2,2,2) + 9*mzv(1,1,4) - mzv(1,5) + 2*mzv(2,1,3))
    print(-mzv(2,4) + 12*mzv(1,1,1,3) - 2*mzv(1,3,2) + 2*mzv(1,1,2,2) - 2*mzv(2,2,2) - 2*mzv(1,1,4))

    -1152*Index(1,7) + 1152*Index(1,1,6) + 9*Index(8)
    4032*Index(1,7) - 2304*Index(1,1,6) - 47*Index(8) + 864*Index(2,6)
    -2304*Index(1,7) - 1728*Index(1,1,6) + 61*Index(8) - 2160*Index(2,6)
    -1728*Index(1,7) + 10368*Index(1,1,6) - 87*Index(8) + 5184*Index(2,6)
    3*Index(8)
    3456*Index(1,7) + 576*Index(1,1,6) - 59*Index(8) + 1728*Index(2,6)
    -5184*Index(1,7) - 3456*Index(1,1,6) + 121*Index(8) - 3888*Index(2,6)
    10368*Index(1,7) + 6912*Index(1,1,6) - 227*Index(8) + 7776*Index(2,6)
    864*Index(1,7) - 5184*Index(1,1,6) + 51*Index(8) - 2592*Index(2,6)
    15*Index(8)
    2016*Index(1,7) + 2304*Index(1,1,6) - 41*Index(8) + 2160*Index(2,6)
    -7200*Index(1,7) - 1728*Index(1,1,6) + 171*Index(8) - 3888*Index(2,6)
    12096*Index(1,7) + 4608*Index(1,1,6) - 153*Index(8) + 7776*Index(2,6)
    -6912*Index(1,7) - 2304*Index(1,1,6) + 455*Index(8) - 4320*Index(2,6)

    [ 3456,576,1728],[-5184,3456,-3888],[ 10368,6912,7776],[ 864,- 5184,- 2592]
=#



@noinline function mzvbasis(n::Int)
    # 1. & 2. インデックスの準備 (変更なし)
    admissible_indexes = admissible_indices(n)
    non_reduction_indexes = Vector{Vector{Int}}()
    reduction_indexes = Vector{Vector{Int}}()

    for idx in admissible_indexes
        if idx ∉ reduction_indexes && collect(dual(IndexWord(idx))) ∉ reduction_indexes
            push!(reduction_indexes, idx)
        else
            push!(non_reduction_indexes, idx)
        end
    end

    sort!(reduction_indexes, by = x->(length(x), x), rev = true)
    admisdict = Dict{IndexWord, Int}(IndexWord(reduction_indexes[i]) => i for i in eachindex(reduction_indexes))
    ncols = length(reduction_indexes)
    
    # 3. & 4. 行列の構築
    raw_rels = relationships_generation(n, admisdict, reduction_indexes)
    nrows = length(raw_rels)
    
    # Bareiss のために BigInt で初期化
    A = zeros(BigInt, nrows, ncols) 
    
    for (r_idx, rel) in enumerate(raw_rels)
        for (word, coeff) in rel.terms
            c_idx = get(admisdict, word, nothing)
            if !isnothing(c_idx)
                # もし coeff が分数なら、ここで整数にスケーリングする必要がある
                # 現状は Integer と想定
                A[r_idx, c_idx] += BigInt(coeff)
            end
        end
    end

    # 5. Bareiss Algorithm (Gauss-Jordan variant)
    pivots = Vector{Int}() 
    curr_row = 1
    prev_pivot = BigInt(1)

    # cntA = 0
    # cntB = 0
    # cntC = 0
    # cntD = 0
    # cntE = 0
    # cntF = 0

    for j in 1:ncols
        # ピボット探索
        pivot_row = findfirst(i -> A[i, j] != 0, curr_row:nrows)
        if isnothing(pivot_row)
            continue 
        end
        pivot_row += (curr_row - 1)
        
        # 行入れ替え
        if curr_row != pivot_row
            A[curr_row, :], A[pivot_row, :] = A[pivot_row, :], A[curr_row, :]
        end
        
        pivot_val = A[curr_row, j]
        push!(pivots, j)
        
        # 他の行の掃き出し
        for i in 1:nrows
            if i == curr_row
                continue
            end

            # A[i, j] が 0 でも、Bareiss の「スケール更新」は必要
            factor = A[i, j] 
            
            for k in 1:ncols

                # cntA += 1
                # (pivot_val      == 0) && (cntB += 1)
                # (A[i,k]         == 0) && (cntC += 1)
                # (factor         == 0) && (cntD += 1)
                # (A[curr_row, k] == 0) && (cntE += 1)
                # (factor == 0) && (A[curr_row, k] == 0) && (cntF += 1)

                num = pivot_val * A[i, k] - factor * A[curr_row, k]
                # ここで必ず割り切れるようになる
                A[i, k] = div(num, prev_pivot) 
            end
            A[i, j] = 0
        end
        
        # 次のステップのために現在のピボットを保存
        prev_pivot = pivot_val
        curr_row += 1
        curr_row > nrows && break
    end

    # println(cntB/cntA,"  ",cntC/cntA,"  ",cntD/cntA,"  ",cntE/cntA,"  ",cntF/cntA,"\n")

    # 5. ループ終了時の prev_pivot が「全行共通の正しい分母」になる
    final_det = prev_pivot 

    # 6. 結果の整理
    basis_indices = setdiff(1:ncols, pivots)
    solutions = Dict{IndexWord, Index}()
    
    for (idx, p_col) in enumerate(pivots)
        # 各行の A[idx, p_col] ではなく、final_det を分母にする
        denom = final_det
        expr = Index()
        for j in (p_col + 1):ncols
            if A[idx, j] != 0
                # A[idx, j] は final_det のスケールまで進化しているので、これで正しい比率になる
                val = Rational{BigInt}(A[idx, j], denom)
                expr -= copy(IndexWord(reduction_indexes[j]) * val)
            end
        end
        # もし Bareiss の過程でピボットの符号が反転している場合、
        # A[idx, p_col] と final_det の符号を比較して調整が必要な場合があるが、
        # 通常はこのまま A[idx, j] / final_det で RREF と一致する。
        solutions[IndexWord(reduction_indexes[p_col])] = expr
    end

    return solutions, [reduction_indexes[b] for b in basis_indices]
end



function relationships_generation(n::Int,admisdict::Dict{IndexWord, Int}, reduction_indexes::Vector{Vector{Int}})
    rels = Vector{Index}()

    # dualでreductionされるものを全てreductionする関数
    function dual_reduction(x::Index)
        r = copy(x)
        for (w,c) in x.terms
            if collect(w) ∉ reduction_indexes
                r -= w*c
                r += dual(w)*c
            end
        end
        return r
    end

    for ad_is in admissible_indices(n-1)
    
        rel = dual_reduction(hoffman_relation(ad_is))
        if rel != 0
            push!(rels,-rel)
        end
    end


    sort!(rels, by = x->(sort([admisdict[idx] for idx in keys(x.terms)])), rev = false)

    deleteat!(rels,2:2:lastindex(rels))


    for w2 in admissible_indices(n-2)
        rel = dual_reduction(f(index(2),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
    end

    for w2 in admissible_indices(n-3)
        rel = dual_reduction(f(index(3),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
        rel = dual_reduction(f(index(1,2),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
    end

    sort!(rels, by = x->(sort([admisdict[idx] for idx in keys(x.terms)])), rev = false)

    return rels
end



@noinline function mzvbasis2(n::Int)
    # 1. & 2. インデックスの準備 (変更なし)
    admissible_indexes = admissible_indices(n)
    non_reduction_indexes = Vector{Vector{Int}}()
    reduction_indexes = Vector{Vector{Int}}()

    for idx in admissible_indexes
        if idx ∉ reduction_indexes && collect(dual(IndexWord(idx))) ∉ reduction_indexes
            push!(reduction_indexes, idx)
        else
            push!(non_reduction_indexes, idx)
        end
    end

    sort!(reduction_indexes, by = x->(length(x), x), rev = true)
    admisdict = Dict{IndexWord, Int}(IndexWord(reduction_indexes[i]) => i for i in eachindex(reduction_indexes))
    ncols = length(reduction_indexes)
    
    # 3. & 4. 行列の構築
    raw_rels = relationships_generation2(n, admisdict, reduction_indexes)
    nrows = length(raw_rels)
    
    # Bareiss のために BigInt で初期化
    #A = zeros(BigInt, nrows, ncols) 
    A = [BigInt() for i in 1:nrows, j in 1:ncols]

    for (r_idx, rel) in enumerate(raw_rels)
        for (word, coeff) in rel.terms
            c_idx = get(admisdict, word, nothing)
            if !isnothing(c_idx)
                # もし coeff が分数なら、ここで整数にスケーリングする必要がある
                # 現状は Integer と想定
                A[r_idx, c_idx] += BigInt(coeff)
            end
        end
    end

    # 5. Bareiss Algorithm (Gauss-Jordan variant)
    pivots = Vector{Int}() 
    curr_row = 1
    prev_pivot = BigInt(1)

    for j in 1:ncols
        # ピボット探索
        pivot_row = findfirst(i -> A[i, j] != 0, curr_row:nrows)
        if isnothing(pivot_row)
            continue 
        end
        pivot_row += (curr_row - 1)
        
        # 行入れ替え
        if curr_row != pivot_row
            A[curr_row, :], A[pivot_row, :] = A[pivot_row, :], A[curr_row, :]
        end
        
        pivot_val = A[curr_row, j]
        push!(pivots, j)
        
        # 他の行の掃き出し
        for i in 1:nrows
            if i == curr_row
                continue
            end
            # factor = A[i, j]
            # if factor == 0
            #     # スケール更新だけは必須（Bareissの整数性を保つため）
            #     for k in (j + 1):ncols
            #         # A[i, k]が0なら計算不要（cntC/cntA 0.68 も効いてくる）
            #         if A[i, k] != 0
            #             A[i, k] = div(pivot_val * A[i, k], prev_pivot)
            #         end
            #     end
            # else
            #     # 通常のBareiss更新
            #     for k in (j + 1):ncols
            #         ak = A[i, k]
            #         ark = A[curr_row, k]
                    
            #         # 統計的に ark == 0 の確率が高いのでここで分岐
            #         if ark == 0
            #             if ak != 0
            #                 A[i, k] = div(pivot_val * ak, prev_pivot)
            #             end
            #         else
            #             # 両方非ゼロ、あるいは ak が 0 のケース
            #             # ここが一番重いパス
            #             A[i, k] = div(pivot_val * ak - factor * ark, prev_pivot)
            #         end
            #     end
            # end

            # A[i, j] が 0 でも、Bareiss の「スケール更新」は必要
            factor = A[i, j]
            if factor == 0
                
                for k in 1:ncols
                    Aik = A[i, k]
                    if iszero(Aik)
                        A[i, k] = BigInt(0)
                    else
                        A[i, k] = div(pivot_val * Aik, prev_pivot) 
                    end
                end
            else
                for k in 1:ncols
                    Aik = A[i,k]
                    Ack = A[curr_row, k]
                    if iszero(Aik)
                        if iszero(Ack)
                            A[i, k] = BigInt(0)
                        else
                            A[i, k] = div(- factor * Ack, prev_pivot)
                        end
                    else
                        if iszero(Ack)
                            A[i, k] = div(pivot_val * Aik, prev_pivot)
                        else
                            A[i, k] = div(pivot_val * Aik - factor * Ack, prev_pivot)
                        end
                    end
                end
            end
            
            A[i, j] = 0
        end
        
        # 次のステップのために現在のピボットを保存
        prev_pivot = pivot_val
        curr_row += 1
        curr_row > nrows && break
    end


    # 5. ループ終了時の prev_pivot が「全行共通の正しい分母」になる
    final_det = prev_pivot 

    # 6. 結果の整理
    basis_indices = setdiff(1:ncols, pivots)
    solutions = Dict{IndexWord, Index}()
    
    for (idx, p_col) in enumerate(pivots)
        # 各行の A[idx, p_col] ではなく、final_det を分母にする
        denom = final_det
        expr = Index()
        for j in (p_col + 1):ncols
            if A[idx, j] != 0
                # A[idx, j] は final_det のスケールまで進化しているので、これで正しい比率になる
                val = Rational{BigInt}(A[idx, j], denom)
                expr -= copy(IndexWord(reduction_indexes[j]) * val)
            end
        end
        # もし Bareiss の過程でピボットの符号が反転している場合、
        # A[idx, p_col] と final_det の符号を比較して調整が必要な場合があるが、
        # 通常はこのまま A[idx, j] / final_det で RREF と一致する。
        solutions[IndexWord(reduction_indexes[p_col])] = expr
    end

    return solutions, [reduction_indexes[b] for b in basis_indices]
end



function relationships_generation2(n::Int,admisdict::Dict{IndexWord, Int}, reduction_indexes::Vector{Vector{Int}})
    rels = Vector{Index}()

    # dualでreductionされるものを全てreductionする関数
    function dual_reduction(x::Index)
        r = copy(x)
        for (w,c) in x.terms
            if collect(w) ∉ reduction_indexes
                r -= w*c
                r += dual(w)*c
            end
        end
        return r
    end

    for ad_is in admissible_indices(n-1)
    
        rel = dual_reduction(hoffman_relation(ad_is))
        if rel != 0
            push!(rels,-rel)
        end
    end


    sort!(rels, by = x->(sort([admisdict[idx] for idx in keys(x.terms)])), rev = false)

    deleteat!(rels,2:2:lastindex(rels))


    for w2 in admissible_indices(n-2)
        rel = dual_reduction(f(index(2),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
    end

    for w2 in admissible_indices(n-3)
        rel = dual_reduction(f(index(3),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
        rel = dual_reduction(f(index(1,2),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
    end

    #sort!(rels, by = x->(sort([admisdict[idx] for idx in keys(x.terms)])), rev = false)

    return rels
end







####################################################################################################################################################################################################









function index_reducing(admissible_indexes)
    
    non_reduction_indexes = Vector{Vector{Int}}()
    reduction_indexes = Vector{Vector{Int}}()

    for idx in admissible_indexes
        if idx ∉ reduction_indexes && collect(dual(IndexWord(idx))) ∉ reduction_indexes
            push!(reduction_indexes, idx)
        else
            push!(non_reduction_indexes, idx)
        end
    end
    
    sort!(reduction_indexes, by = x->(length(x), x), rev = true)

    return reduction_indexes, non_reduction_indexes
end


function relationships_generation3(n::Int,admisdict::Dict{IndexWord, Int}, reduction_indexes::Vector{Vector{Int}})
    rels = Vector{Index}()

    # dualでreductionされるものを全てreductionする関数
    function dual_reduction(x::Index)
        r = copy(x)
        for (w,c) in x.terms
            if collect(w) ∉ reduction_indexes
                r -= w*c
                r += dual(w)*c
            end
        end
        return r
    end

    for ad_is in admissible_indices(n-1)
    
        rel = dual_reduction(hoffman_relation(ad_is))
        if rel != 0
            push!(rels,-rel)
        end
    end


    sort!(rels, by = x->(sort([admisdict[idx] for idx in keys(x.terms)])), rev = false)

    deleteat!(rels,2:2:lastindex(rels))


    for w2 in admissible_indices(n-2)
        rel = dual_reduction(f(index(2),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
    end

    for w2 in admissible_indices(n-3)
        rel = dual_reduction(f(index(3),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
        rel = dual_reduction(f(index(1,2),index(w2)))
        if rel != 0
            push!(rels,rel)
        end
    end

    #sort!(rels, by = x->(sort([admisdict[idx] for idx in keys(x.terms)])), rev = false)

    return rels
end


function build_matrix(raw_rels::Vector{Index}, admisdict, nrows, ncols)
    # Bareiss のために BigInt で初期化 
    A = [BigInt() for i in 1:nrows, j in 1:ncols]

    #
    for (r_idx, rel) in enumerate(raw_rels)
        for (word, coeff) in rel.terms
            c_idx = get(admisdict, word, nothing)
            if !isnothing(c_idx)
                # coeffは整数
                A[r_idx, c_idx] += BigInt(coeff)
            end
        end
    end
    return A
end

function bareiss_algorithm(A, nrows, ncols)
    # 5. Bareiss Algorithm (Gauss-Jordan variant)
    pivots = Vector{Int}() 
    curr_row = 1
    prev_pivot = BigInt(1)

    for j in 1:ncols
        # ピボット探索
        pivot_row = findfirst(i -> A[i, j] != 0, curr_row:nrows)
        if isnothing(pivot_row)
            continue 
        end
        pivot_row += (curr_row - 1)
        
        # 行入れ替え
        if curr_row != pivot_row
            A[curr_row, :], A[pivot_row, :] = A[pivot_row, :], A[curr_row, :]
        end
        
        pivot_val = A[curr_row, j]
        push!(pivots, j)
        
        # 他の行の掃き出し
        for i in 1:nrows
            if i == curr_row
                continue
            end

            # A[i, j] が 0 でも、Bareiss の「スケール更新」は必要
            factor = A[i, j]
            if factor == 0
                
                for k in 1:ncols
                    Aik = A[i, k]
                    if iszero(Aik)
                        A[i, k] = BigInt(0)
                    else
                        A[i, k] = div(pivot_val * Aik, prev_pivot) 
                    end
                end
            else
                for k in 1:ncols
                    Aik = A[i,k]
                    Ack = A[curr_row, k]
                    if iszero(Aik)
                        if iszero(Ack)
                            A[i, k] = BigInt(0)
                        else
                            A[i, k] = div(- factor * Ack, prev_pivot)
                        end
                    else
                        if iszero(Ack)
                            A[i, k] = div(pivot_val * Aik, prev_pivot)
                        else
                            A[i, k] = div(pivot_val * Aik - factor * Ack, prev_pivot)
                        end
                    end
                end
            end
            
            A[i, j] = 0
        end
        
        # 次のステップのために現在のピボットを保存
        prev_pivot = pivot_val
        curr_row += 1
        curr_row > nrows && break
    end
    return pivots, prev_pivot
end

@noinline function mzvbasis3(n::Int)

    admissible_indexes = admissible_indices(n)


    reduction_indexes, non_reduction_indexes = index_reducing(admissible_indexes)


    # 列番号 => 具体的なインデックス の辞書
    admisdict = Dict{IndexWord, Int}(IndexWord(reduction_indexes[i]) => i for i in eachindex(reduction_indexes))
    # メイン行列の列数
    ncols = length(reduction_indexes)
    
    
    # raw_rels::Vector{Index} で、raw_rels[i] = 0 となるIndex関係式の集まり
    raw_rels = relationships_generation3(n, admisdict, reduction_indexes)
    # メイン行列の行数(関係式の数)
    nrows = length(raw_rels)
    
    
    A = build_matrix(raw_rels, admisdict, nrows, ncols)

    pivots, final_det = bareiss_algorithm(A, nrows, ncols)


    # 6. 結果の整理
    basis_indices = setdiff(1:ncols, pivots)
    solutions = Dict{IndexWord, Index}()
    
    for (idx, p_col) in enumerate(pivots)
        # 各行の A[idx, p_col] ではなく、final_det を分母にする
        denom = final_det
        expr = Index()
        for j in (p_col + 1):ncols
            if A[idx, j] != 0
                # A[idx, j] は final_det のスケールまで進化しているので、これで正しい比率になる
                val = Rational{BigInt}(A[idx, j], denom)
                expr -= copy(IndexWord(reduction_indexes[j]) * val)
            end
        end
        # もし Bareiss の過程でピボットの符号が反転している場合、
        # A[idx, p_col] と final_det の符号を比較して調整が必要な場合があるが、
        # 通常はこのまま A[idx, j] / final_det で RREF と一致する。
        solutions[IndexWord(reduction_indexes[p_col])] = expr
    end

    return solutions, [reduction_indexes[b] for b in basis_indices]
end