#[ monomials.jl ]#

# This file defines functions related to monomials

#=
export shuffle_product, stuffle_product, star_stuffle_product, 
       shuffle_product_double, stuffle_product_double, star_stuffle_product_double,
       shuffle_pow, stuffle_pow, star_stuffle_pow, 
       shpw, stpw, starstpw,
       st_index1_pow, sh_index1_pow,
       ⩊, ∗, ⋆,
       Hoffman_hom, Hoffman_antihom, starword_to_word,
       dual, Hoffman_dual, Landen_dual,
       stuffle_regularization_polynomial, shuffle_regularization_polynomial,
       rho_t, rho, reg_st, reg_sh,
=#

"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""

###################################################################################################
############## monomial functions #################################################################

#=
last_letter(w::Word)::Word = length(w) == 0 ? () : Word(w[end])
remove_last_letter(w::Word)::Word = w[1:end-1]
last_lettering(w::Word, g::ExprInt)::Int64 = (ig = findlast(==(g),w)) |> isnothing ? 0 : ig
evaporate_last_lettering(w::Word, g::ExprInt)::Word = w[last_lettering(w,g):end]
remove_last_letter(w::Word, g::ExprInt)::Word = w[1:last_lettering(w,g)]

function last_letter_idx(v::Word)::Word
    if length(v) == 0
        return ()
    elseif v[end] == 1
        return (2)
    else
        return (1)
    end
end
function remove_last_letter_idx(v::Word)::Word
    if v[end] == 1
        return v[1:end-1]
    else
        r = copy(v)
        r[end] -= 1
        return r
    end
end
function last_letter_v(v::Vector{Int})::Vector{Int}
    if lastindex(v) == 0
        return Int[]
    else
        return [v[end]]
    end
end
function remove_last_letter_v(v::Vector{Int})::Vector{Int}
    return v[1:end-1]
end
=#

#last_lettering(v::Idv, g::ExprInt)::Idv = (ig = findlast(==(g),v)) |> isnothing ? 0 : ig

###################################################################################################
############## monomial products ##################################################################

function monomial_sh(a::Vector{Int},b::Vector{Int})::Vector{Vector{Int}}
    lb = lastindex(b)
    la = lastindex(a)
    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[1:i]]
    end

    for i in 1:lb
        
        h1[1] = [b[1:i]]

        for j in 1:la
            h1[j+1] = vcat(vcat.(h1[j],a[j]),vcat.(h1[j+1],b[i]))
        end
    end

    return h1[end]
end
function monomial_sh_r(a::Vector{Int},b::Vector{Int})::Vector{Vector{Int}}
    lb = lastindex(b)
    la = lastindex(a)
    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[la+1-i:la]]
    end

    for i in 1:lb
        
        h1[1] = [b[lb+1-i:lb]]

        for j in 1:la
            h1[j+1] = vcat(vcat.(a[la+1-j],h1[j]),vcat.(b[lb+1-i],h1[j+1]))
        end
    end

    return h1[end]
end
function monomial_st(a::Vector{Int},b::Vector{Int})::Vector{Vector{Int}}
    lb = lastindex(b)
    la = lastindex(a)
    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[1:i]]
    end

    for i in 1:lb

        buf1 = [b[1:i]]

        for j in 1:la
            buf2 = h1[j]
            h1[j] = buf1
            buf1 = vcat(vcat.(buf1,a[j]),vcat.(h1[j+1],b[i]),vcat.(buf2,a[j]+b[i]))
        end

        h1[end] = buf1
    end

    return h1[end]
end
function monomial_st_r(a::Vector{Int},b::Vector{Int})::Vector{Vector{Int}}
    lb = lastindex(b)
    la = lastindex(a)
    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[la+1-i:la]]
    end

    for i in 1:lb

        buf1 = [b[lb+1-i:lb]]

        for j in 1:la
            buf2 = h1[j]
            h1[j] = buf1
            buf1 = vcat(vcat.(a[la+1-j],buf1),vcat.(b[lb+1-i],h1[j+1]),vcat.(a[la+1-j]+b[lb+1-i],buf2))
        end

        h1[end] = buf1
    end

    return h1[end]
end
function monomial_st_star(a::Vector{Int}, b::Vector{Int})::Tuple{Vector{Vector{Int}}, Vector{Bool}}
    # 通常の stuffle（符号なし）
    words = monomial_st(a, b)

    dep_sum = (lastindex(a) + lastindex(b)) & 1
    n = lastindex(words)
    signs = Vector{Bool}(undef, n)

    for i in 1:n
        dep = lastindex(words[i])
        # 偶奇が同じなら +（false）、異なれば −（true）
        signs[i] = xor(dep_sum,dep) & 1
    end

    return words, signs
end
function monomial_st_star_r(a::Vector{Int}, b::Vector{Int})::Tuple{Vector{Vector{Int}}, Vector{Bool}}
    # 通常の stuffle（符号なし）
    words = monomial_st_r(a, b)

    dep_sum = (lastindex(a) + lastindex(b)) & 1
    n = lastindex(words)
    signs = Vector{Bool}(undef, n)

    for i in 1:n
        dep = lastindex(words[i])
        # 偶奇が同じなら +（false）、異なれば −（true）
        signs[i] = xor(dep_sum,dep) & 1
    end

    return words, signs
end

function monomial_sh_double(a::Vector{Int})::Vector{Vector{Int}}
    la = lastindex(a)

    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[1:i]]
    end

    for i in 1:la
        h1[i] = h1[i+1]

        for j in i:la
            h1[j+1] = vcat(vcat.(h1[j],a[j]),vcat.(h1[j+1],a[i]))
        end
    end

    return h1[end]
end
function monomial_sh_double_r(a::Vector{Int})::Vector{Vector{Int}}
    la = lastindex(a)

    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[la+1-i:la]]
    end

    for i in 1:la
        h1[i] = h1[i+1]

        for j in i:la
            h1[j+1] = vcat(vcat.(a[la+1-j],h1[j]),vcat.(a[la+1-i],h1[j+1]))
        end
    end

    return h1[end]
end
function monomial_st_double(a::Vector{Int})::Vector{Vector{Int}}
    la = lastindex(a)
    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[1:i]]
    end

    for i in 1:la

        buf1 = h1[i+1]

        for j in i:la
            buf2 = h1[j]
            h1[j] = buf1
            buf1 = vcat(vcat.(buf1,a[j]),vcat.(h1[j+1],a[i]),vcat.(buf2,a[j]+a[i]))
        end

        h1[end] = buf1
    end

    return h1[end]
end
function monomial_st_double_r(a::Vector{Int})::Vector{Vector{Int}}
    la = lastindex(a)
    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[la+1-i:la]]
    end

    for i in 1:la

        buf1 = h1[i+1]

        for j in i:la
            buf2 = h1[j]
            h1[j] = buf1
            buf1 = vcat(vcat.(a[la+1-j],buf1),vcat.(a[la+1-i],h1[j+1]),vcat.(a[la+1-j]+a[la+1-i],buf2))
        end

        h1[end] = buf1
    end

    return h1[end]
end
function monomial_st_star_double(a::Vector{Int})::Tuple{Vector{Vector{Int}},Vector{Bool}}
    # 通常の stuffle（符号なし）
    words = monomial_st_double(a)

    n = lastindex(words)
    signs = Vector{Bool}(undef, n)
    
    for i in 1:n
        # 偶奇が同じなら +（false）、異なれば −（true）
        signs[i] = lastindex(words[i]) & 1
    end
    
    return words, signs
end
function monomial_st_star_double_r(a::Vector{Int})::Tuple{Vector{Vector{Int}},Vector{Bool}}
    # 通常の stuffle（符号なし）
    words = monomial_st_double_r(a)

    n = lastindex(words)
    signs = Vector{Bool}(undef, n)
    
    for i in 1:n
        # 偶奇が同じなら +（false）、異なれば −（true）
        signs[i] = lastindex(words[i]) & 1
    end
    
    return words, signs
end


###################################################################################################
############## products ###########################################################################

function shuffle_product(a::Hoffman, b::Hoffman)::Hoffman
    s1 = Hoffman()
    if get_index_orientation()
        for (ma,ca) in a.terms
            s2 = Hoffman()
            for (mb,cb) in b.terms
                add!(s2,Hoffman(monomial_sh(ma.tovec,mb.tovec)),cb)
            end
            add!(s1,s2,ca)
        end
    else
        for (ma,ca) in a.terms
            s2 = Hoffman()
            for (mb,cb) in b.terms
                add!(s2,Hoffman(monomial_sh_r(ma.tovec,mb.tovec)),cb)
            end
            add!(s1,s2,ca)
        end
    end

    return s1
end
function stuffle_product(a::Index, b::Index)::Index
    s1 = Index()
    if get_index_orientation()
        for (ma,ca) in a.terms
            s2 = Index()
            for (mb,cb) in b.terms
                add!(s2,Index(monomial_st(ma.tovec,mb.tovec)),cb)
            end
            add!(s1,s2,ca)
        end
    else
        for (ma,ca) in a.terms
            s2 = Index()
            for (mb,cb) in b.terms
                add!(s2,Index(monomial_st_r(ma.tovec,mb.tovec)),cb)
            end
            add!(s1,s2,ca)
        end
    end
    return s1
end
function star_stuffle_product(a::Index, b::Index)::Index
    s1 = Index()
    if get_index_orientation()
        for (ma,ca) in a.terms
            s2 = Index()
            for (mb,cb) in b.terms
                add!(s2,Index(monomial_st_star(ma.tovec,mb.tovec)),cb)
            end
            add!(s1,s2,ca)
        end
    else
        for (ma,ca) in a.terms
            s2 = Index()
            for (mb,cb) in b.terms
                add!(s2,Index(monomial_st_star_r(ma.tovec,mb.tovec)),cb)
            end
            add!(s1,s2,ca)
        end
    end
    return s1
end
function shuffle_product_double(a::Hoffman)::Hoffman
    result = Hoffman()
    pairs = collect(a.terms)
    n = lastindex(pairs)

    if get_index_orientation()
        for i in 1:n
            wi, ci = pairs[i]
            vi = collect(wi)
    
            # 自乗項
            add!(result, Hoffman(monomial_sh_double(vi)), ci^2)
    
            # 交差項
            for j in i+1:n
                wj, cj = pairs[j]
                vj = collect(wj)
    
                add!(result, Hoffman(monomial_sh(vi, vj)), 2*ci*cj)
            end
    
        end
    else
        for i in 1:n
            wi, ci = pairs[i]
            vi = collect(wi)
    
            # 自乗項
            add!(result, Hoffman(monomial_sh_double_r(vi)), ci^2)
    
            # 交差項
            for j in i+1:n
                wj, cj = pairs[j]
                vj = collect(wj)
    
                add!(result, Hoffman(monomial_sh_r(vi, vj)), 2*ci*cj)
            end
    
        end
    end

    return result
end
function stuffle_product_double(a::Index)::Index
    result = Index()
    pairs = collect(a.terms)
    n = lastindex(pairs)

    if get_index_orientation()
        for i in 1:n
            wi, ci = pairs[i]
            vi = collect(wi)
    
            # 自乗項
            add!(result, Index(monomial_st_double(vi)), ci^2)
    
            # 交差項
            for j in i+1:n
                wj, cj = pairs[j]
                vj = collect(wj)
    
                add!(result, Index(monomial_st(vi, vj)), 2*ci*cj)
            end
        end
    else
        for i in 1:n
            wi, ci = pairs[i]
            vi = collect(wi)
    
            # 自乗項
            add!(result, Index(monomial_st_double_r(vi)), ci^2)
    
            # 交差項
            for j in i+1:n
                wj, cj = pairs[j]
                vj = collect(wj)
    
                add!(result, Index(monomial_st_r(vi, vj)), 2*ci*cj)
            end
        end
    end

    return result
end
function star_stuffle_product_double(a::Index)::Index
    result = Index()
    pairs = collect(a.terms)
    n = lastindex(pairs)

    if get_index_orientation()
        for i in 1:n
            wi, ci = pairs[i]
            vi = collect(wi)
    
            # 自乗項
            add!(result, Index(monomial_st_star_double(vi)), ci^2)
    
            # 交差項
            for j in i+1:n
                wj, cj = pairs[j]
                vj = collect(wj)
    
                add!(result, Index(monomial_st_star(vi, vj)), 2*ci*cj)
            end
        end
    else
        for i in 1:n
            wi, ci = pairs[i]
            vi = collect(wi)
    
            # 自乗項
            add!(result, Index(monomial_st_star_double_r(vi)), ci^2)
    
            # 交差項
            for j in i+1:n
                wj, cj = pairs[j]
                vj = collect(wj)
    
                add!(result, Index(monomial_st_star_r(vi, vj)), 2*ci*cj)
            end
        end
    end

    return result
end

# 一般的に早い
function shuffle_pow(a::Hoffman,n::Int)::Hoffman
    base = deepcopy(a)
    r = one(Hoffman)

    if n < 0
        throw(ArgumentError("shuffle_pow is only defined for nonnegative powers"))
    elseif n == 0
        return r
    elseif n == 1
        return base
    end

    if n&1 == 1
        r = shuffle_product(r,base)
    end
    n >>= 1
    while n>0
        base = shuffle_product_double(base)
        if n&1 == 1
            r = shuffle_product(r,base)
        end
        n>>=1
    end

    return r
end
function stuffle_pow(a::Index,n::Int)::Index
    base = deepcopy(a)
    r = one(Index)

    if n < 0
        throw(ArgumentError("stuffle_pow is only defined for nonnegative powers"))
    elseif n == 0
        return r
    elseif n == 1
        return base
    end
    
    if n&1 == 1
        r = stuffle_product(r,base)
    end
    n >>= 1
    while n>0
        base = stuffle_product_double(base)
        if n&1 == 1
            r = stuffle_product(r,base)
        end
        n>>=1
    end
    return r
end
function star_stuffle_pow(a::Index,n::Int)::Index
    base = deepcopy(a)
    r = one(Index)

    if n < 0
        throw(ArgumentError("stuffle_pow is only defined for nonnegative powers"))
    elseif n == 0
        return r
    elseif n == 1
        return base
    end
    
    if n&1 == 1
        r = star_stuffle_product(r,base)
    end
    n >>= 1
    while n>0
        base = star_stuffle_product_double(base)
        if n&1 == 1
            r = star_stuffle_product(r,base)
        end
        n>>=1
    end
    return r
end
# aの項が2項とかの時に速い
function shpw(a::Hoffman,n::Int)::Hoffman
    r = one(Hoffman)
    if n==0
        return r
    else
        return shuffle_product(shpw(a,n-1),a)
    end
end
function stpw(a::Index,n::Int)::Index
    r = one(Index)
    if n==0
        return r
    else
        return stuffle_product(stpw(a,n-1),a)
    end
end
function starstpw(a::Index,n::Int)::Index
    r = one(Index)
    if n==0
        return r
    else
        return star_stuffle_product(starstpw(a,n-1),a)
    end
end

"""
    st_index1_pow(n::Int) -> Index

st_index1_pow(n) == stpw(Index(1),n)
"""
function st_index1_pow(n::Int)::Index
    if n == 0
        return one(Index)
    end
    # 返す値
    idx = Index()

    # 0 を許さず、最上位ビットだけ 1 に固定
    # 残り n-1 個のビットを全探索
    total = UInt(1) << (n-1)

    # n! を BigInt で計算
    nf = factorial(BigInt(n))

    for mask in UInt(0):(total-1)
        # ビット列を作成：必ず先頭は 1
        bits = ones(Int,n)
        bits[1] = 2

        # 残りのビット（下位 n-1 を mask から読む）
        for i in 2:n
            bits[i] += (mask >> (n-i)) & 0x1 == 1 ? 1 : 0
        end

        # 100101 → [3,2,1]
        v = idxprs(Word(bits))
        
        # multinomial
        coeff = multinomial(v,nf)

        #@show (bits,v,w,coeff)

        idx.terms[Word(v)] = coeff
    end

    return idx
end

function sh_y_pow(n::Int)::Hoffman
    h = Hoffman()
    h.terms[Word(fill(2,n))] = factorial(BigInt(n))
    return h
end

shuffle_product(a::Index, b::Index)::Index            = Index(shuffle_product(a.toHoffman, b.toHoffman))
shuffle_product_double(a::Index)::Index               = Index(shuffle_product_double(a.toHoffman))
shuffle_pow(a::Index, n::Int)::Index                  = Index(shuffle_pow(a.toHoffman, n))
shpw(a::Index, n::Int)::Index                         = Index(shpw(a.toHoffman, n))

stuffle_product(a::Hoffman, b::Hoffman)::Hoffman      = Hoffman(stuffle_product(a.toIndex, b.toIndex))
stuffle_product_double(a::Hoffman)::Hoffman           = Hoffman(stuffle_product_double(a.toIndex))
stuffle_pow(a::Hoffman, n::Int)::Hoffman              = Hoffman(stuffle_pow(a.toIndex, n))
stpw(a::Hoffman, n::Int)::Hoffman                     = Hoffman(stpw(a.toIndex, n))

star_stuffle_product(a::Hoffman, b::Hoffman)::Hoffman = Hoffman(star_stuffle_product(a.toIndex, b.toIndex))
star_stuffle_product_double(a::Hoffman)::Hoffman      = Hoffman(star_stuffle_product_double(a.toIndex))
star_stuffle_pow(a::Hoffman, n::Int)::Hoffman         = Hoffman(star_stuffle_pow(a.toIndex, n))
starstpw(a::Hoffman, n::Int)::Hoffman                 = Hoffman(starstpw(a.toIndex, n))

⩊(a,b) = shuffle_product(a,b)
∗(a,b) = stuffle_product(a,b)
⋆(a,b) = star_stuffle_product(a,b)

###################################################################################################
############## homomorphic ########################################################################

@inline function Hoffman_hom(w::Word, image::Vector{Hoffman})::Hoffman
    s = one(Hoffman)
    for idx in w
        s = s * image[idx]
    end
    return s
end
@inline function Hoffman_antihom(w::Word, image::Vector{Hoffman})::Hoffman
    s = one(Hoffman)
    for idx in w
        s = image[idx] * s
    end
    return s
end

#star-valueのindexに対応するwordから通常のindexに
function monomial_sw_w(w::Word)::Hoffman
    s = one(Hoffman)
    for idx in w[2:end]
        if idx == 1
            s *= x
        else
            s *= x+y
        end
    end
    return y * s
end
function monomial_sw_w_r(w::Word)::Hoffman
    s = one(Hoffman)
    for i in lastindex(w)-1:-1:1
        idx = w[i]
        if idx == 1
            s = x * s
        else
            s = (x+y) * s
        end
    end
    return s * y
end

function starword_to_word(w::Hoffman)::Hoffman
    s = Hoffman()
    if get_index_orientation()
        for (mw,cw) in w.terms
            add!(s,monomial_sw_w(mw),cw)
        end
    else
        for (mw,cw) in w.terms
            add!(s,monomial_sw_w_r(mw),cw)
        end
    end

    return s
end
starword_to_word(i::Index)::Index = Index(starword_to_word(i.toHoffman))


###################################################################################################
############## dual index #########################################################################


# 双対インデックス
function monomial_dual_h(w::Word)::Word
    r = Word(Tuple((w[i]==1 ? 2 : 1) for i in lastindex(w):-1:1))
    return r
end
function monomial_dual_i(w::Word)::Word
    lw = lastindex(w)
    l = sum(w) - lw
    v = ones(Int,l)
    cid = 0
    for i in lw:-1:1
        cid += w[i]-1
        v[cid] += 1
    end
    return Word(v)
end
function monomial_dual_i_r(w::Word)::Word
    lw = lastindex(w)
    l = sum(w) - lw
    v = ones(Int,l)
    cid = l+1
    for i in 1:lw
        cid -= w[i]-1
        v[cid] += 1
    end
    return Word(v)
end
function dual(w::Hoffman)::Hoffman
    s = Hoffman()
    for (mw,cw) in w.terms
        mwd = monomial_dual_h(mw)
        if haskey(s.terms,mwd)
            if s.terms[mwd] == -cw
                delete!(s.terms,mwd)
            else
                s.terms[mwd] += cw
            end
        else
            s.terms[mwd] = cw
        end
    end
    return s
end
function dual(i::Index)::Index
    s = Index()
    if get_index_orientation()
        for (mw,cw) in i.terms
            mwd = monomial_dual_i(mw)
            if haskey(s.terms,mwd)
                if s.terms[mwd] == -cw
                    delete!(s.terms,mwd)
                else
                    s.terms[mwd] += cw
                end
            else
                s.terms[mwd] = cw
            end
        end
    else
        for (mw,cw) in i.terms
            mwd = monomial_dual_i_r(mw)
            if haskey(s.terms,mwd)
                if s.terms[mwd] == -cw
                    delete!(s.terms,mwd)
                else
                    s.terms[mwd] += cw
                end
            else
                s.terms[mwd] = cw
            end
        end
    end

    return s
end

# Hoffman双対インデックス
function monomial_hof_dual_h(w::Word)::Word
    lw = lastindex(w)
    if lw == 0
        return Word()
    end
    return Word(2) * Word(Tuple(3-i for i in w[2:lw]))
end
function monomial_hof_dual_h_r(w::Word)::Word
    lw = lastindex(w)
    if lw == 0
        return Word()
    end
    return Word(Tuple(3-i for i in w[1:lw-1])) * Word(2)
end
function Hoffman_dual(w::Hoffman)::Hoffman
    s = Hoffman()
    if get_index_orientation()
        for (mw,cw) in w.terms
            mwd = monomial_hof_dual_h(mw)
            if haskey(s.terms,mwd)
                if s.terms[mwd] == -cw
                    delete!(s.terms,mwd)
                else
                    s.terms[mwd] += cw
                end
            else
                s.terms[mwd] = cw
            end
        end
    else
        for (mw,cw) in w.terms
            mwd = monomial_hof_dual_h_r(mw)
            if haskey(s.terms,mwd)
                if s.terms[mwd] == -cw
                    delete!(s.terms,mwd)
                else
                    s.terms[mwd] += cw
                end
            else
                s.terms[mwd] = cw
            end
        end
    end

    return s
end
function monomial_hof_dual_i(w::Word)::Word
    lw = lastindex(w)
    l = sum(w) - lw + 1
    v = ones(Int,l)
    cid = 1
    for i in 1:lw-1
        cid += w[i]-1
        v[cid] += 1
    end
    return Word(v)
end
# 実は等価
function monomial_hof_dual_i_r(w::Word)::Word
    lw = lastindex(w)
    l = sum(w) - lw + 1
    v = ones(Int,l)
    cid = 1
    for i in 1:lw-1
        cid += w[i]-1
        v[cid] += 1
    end
    return Word(v)
end
# monomial_hof_dual_i と monomial_hof_dual_i_r が等価なので↓は場合分けしなくてよい
function Hoffman_dual(idx::Index)::Index
    s = Index()
    for (mw,cw) in idx.terms
        mwd = monomial_hof_dual_i(mw)
        if haskey(s.terms,mwd)
            if s.terms[mwd] == -cw
                delete!(s.terms,mwd)
            else
                s.terms[mwd] += cw
            end
        else
            s.terms[mwd] = cw
        end
    end
    return s
end

# Landen dual
function monomial_Landen(w::Word)::Hoffman
    s = one(Hoffman)
    c1 = y+x
    for i in 1:lastindex(w)
        if w[i] == 1
            s = s*c1
        else
            s = -(s*y)
        end
    end
    return s
end
function monomial_Landen_r(w::Word)::Hoffman
    s = one(Hoffman)
    c1 = y+x
    for i in lastindex(w):-1:1
        if w[i] == 1
            s = c1*s
        else
            s = -(y*s)
        end
    end
    return s
end
function Landen_dual(w::Hoffman)::Hoffman
    s = Hoffman()
    if get_index_orientation()
        for (mw,cw) in w.terms
            add!(s,monomial_Landen(mw),cw)
        end
    else
        for (mw,cw) in w.terms
            add!(s,monomial_Landen_r(mw),cw)
        end
    end
    return s
end
Landen_dual(idx::Word)::Word = Index(monomial_Landen(IndexWordtoHoffmanWord(idx)))
Landen_dual(idx::Index)::Index = Index(Landen_dual(idx.toHoffman))

###################################################################################################
############## regularization polynomial ##########################################################


# 右から1が最も多く続くIndexを返す
# 同数なら全てを返す
function rightmostones(i::Index)::Vector{MonoIndex}
    maxlen = -1
    winners = Vector{Tuple{Word,Rational{BigInt}}}()

    # 辞書を一度だけ走査
    for (w, c) in i.terms
        # 末尾の1の連続数をカウント
        cnt = 0
        for x in lastindex(w.t):-1:1
            w[x] == 1 ? (cnt += 1) : break
        end

        if cnt > maxlen
            maxlen = cnt
            empty!(winners)
            push!(winners, (w, c))
        elseif cnt == maxlen
            push!(winners, (w, c))
        end
    end

    # 結果をMonoIndex型の配列に
    res = Vector{MonoIndex}(undef, lastindex(winners))
    for (k, (w, c)) in enumerate(winners)
        res[k] = MonoIndex(w,c)
    end
    return res
end
# 右からyが最も多く続くHoffmanを返す
# 同数なら全てを返す
function rightmostys(h::Hoffman)::Vector{MonoIndex}
    maxlen = -1
    winners = Vector{Tuple{Word,Rational{BigInt}}}()

    # 辞書を一度だけ走査
    for (w, c) in h.terms
        # 末尾の2の連続数をカウント
        cnt = 0
        for x in lastindex(w.t):-1:1
            w[x] == 2 ? (cnt += 1) : break
        end

        if cnt > maxlen
            maxlen = cnt
            empty!(winners)
            push!(winners, (w, c))
        elseif cnt == maxlen
            push!(winners, (w, c))
        end
    end

    # 結果をMonoIndex型の配列に
    # MonoIndexをただWordと係数を一度に入れれる便利なTupleのように扱っている(Wordの扱いはHoffmanとしている)
    res = Vector{MonoIndex}(undef, lastindex(winners))
    for (k, (w, c)) in enumerate(winners)
        res[k] = MonoIndex(w,c)
    end
    return res
end
function leftmostones(i::Index)::Vector{MonoIndex}
    maxlen = -1
    winners = Vector{Tuple{Word,Rational{BigInt}}}()

    # 辞書を一度だけ走査
    for (w, c) in i.terms
        # 末尾の1の連続数をカウント
        cnt = 0
        for x in 1:lastindex(w.t)
            w[x] == 1 ? (cnt += 1) : break
        end

        if cnt > maxlen
            maxlen = cnt
            empty!(winners)
            push!(winners, (w, c))
        elseif cnt == maxlen
            push!(winners, (w, c))
        end
    end

    # 結果をMonoIndex型の配列に
    res = Vector{MonoIndex}(undef, lastindex(winners))
    for (k, (w, c)) in enumerate(winners)
        res[k] = MonoIndex(w,c)
    end
    return res
end
function leftmostys(h::Hoffman)::Vector{MonoIndex}
    maxlen = -1
    winners = Vector{Tuple{Word,Rational{BigInt}}}()

    # 辞書を一度だけ走査
    for (w, c) in h.terms
        # 末尾の2の連続数をカウント
        cnt = 0
        for x in 1:lastindex(w.t)
            w[x] == 2 ? (cnt += 1) : break
        end

        if cnt > maxlen
            maxlen = cnt
            empty!(winners)
            push!(winners, (w, c))
        elseif cnt == maxlen
            push!(winners, (w, c))
        end
    end

    # 結果をMonoIndex型の配列に
    # MonoIndexをただWordと係数を一度に入れれる便利なTupleのように扱っている(Wordの扱いはHoffmanとしている)
    res = Vector{MonoIndex}(undef, lastindex(winners))
    for (k, (w, c)) in enumerate(winners)
        res[k] = MonoIndex(w,c)
    end
    return res
end
# 右から続く1の長さとそれ以外の部分をTupleにして返す
function rightonesplit(w::Word)::Tuple{Word,Int}
    cnt = 0
    for x in lastindex(w):-1:1
        w[x] == 1 ? (cnt += 1) : break
    end
    return (w[1:end-cnt] , cnt)
end
function rightysplit(w::Word)::Tuple{Word,Int}
    cnt = 0
    for x in lastindex(w):-1:1
        w[x] == 2 ? (cnt += 1) : break
    end
    return (w[1:end-cnt] , cnt)
end
function leftonesplit(w::Word)::Tuple{Word,Int}
    cnt = 0
    for x in 1:lastindex(w)
        w[x] == 1 ? (cnt += 1) : break
    end
    return (w[cnt+1:end] , cnt)
end
function leftysplit(w::Word)::Tuple{Word,Int}
    cnt = 0
    for x in 1:lastindex(w)
        w[x] == 2 ? (cnt += 1) : break
    end
    return (w[cnt+1:end] , cnt)
end
#=
1. rightmostonesをする
2. 各項に対してshuffle_product(後ろ,1,1,...)をする
=#
function stuffle_regularization_polynomial(w::Word)::Poly{Index}
    isone(w) && return one(Poly{Index})
    if get_index_orientation()
        if w[end] != 1  # xに相当する
            return Poly(IndexWordtoIndex(w))
        end
        z = IndexWordtoIndex(w)
        r = Poly{Index}()

        while !iszero(z)    # z が0でない限り続ける
            m = rightmostones(z)
            zadd = Index()
            mi_d = 0
            for mi in m
                # zのうち右から続く1が一番多いものの集まりm
                # mの一つ一つの1の続く部分(mi_d個)とそれ以外(mi_r)
                mi_r, mi_d = rightonesplit(mi.word)

                # mi_r ∗ Index(1)^{∗mi_d} を計算
                #mi_t = stuffle_product(IndexWordtoIndex(mi_r),stpw(Index(1),mi_d))
                mi_t = stuffle_product(IndexWordtoIndex(mi_r),st_index1_pow(mi_d))

                # 係数を合わせる
                diffc = mi.coeff//mi_t.terms[mi.word]
                
                # zからひとつづつ引く
                z -= mi_t * diffc

                # r(正規化多項式)に足しこむ分
                zadd.terms[mi_r] = diffc
            end
            r.terms[mi_d] = zadd
        end
        return r
    else
        if w[1] != 1
            return Poly(IndexWordtoIndex(w))
        end
        z = IndexWordtoIndex(w)
        r = Poly{Index}()
        
        while !iszero(z)
            m = leftmostones(z)
            zadd = Index()
            mi_d = 0
            for mi in m
                mi_r, mi_d = leftonesplit(mi.word)

                #mi_t = stuffle_product(IndexWordtoIndex(mi_r),stpw(Index(1),mi_d))
                mi_t = stuffle_product(IndexWordtoIndex(mi_r),st_index1_pow(mi_d))
                
                diffc = mi.coeff//mi_t.terms[mi.word]

                z -= mi_t * diffc

                zadd.terms[mi_r] = diffc
            end
            r.terms[mi_d] = zadd
        end
        return r
    end
end
function stuffle_regularization_polynomial(idx::Index)::Poly{Index}
    rp = Poly{Index}()
    for (w,c) in idx.terms
        rp += c * stuffle_regularization_polynomial(w)
    end
    return rp
end

function shuffle_regularization_polynomial(w::Word)::Poly{Hoffman}
    isone(w) && return one(Poly{Hoffman})
    if get_index_orientation()
        if w[end] != 2  # xに相当する
            return Poly(HoffmanWordtoHoffman(w))
        end
        z = HoffmanWordtoHoffman(w)
        r = Poly{Hoffman}()

        while !iszero(z)    # z が0でない限り続ける
            m = rightmostys(z)
            zadd = Hoffman()
            mi_d = 0
            for mi in m
                # zのうち右から続く1が一番多いものの集まりm
                # mの一つ一つの1の続く部分(mi_d個)とそれ以外(mi_r)
                mi_r, mi_d = rightysplit(mi.word)

                # mi_r ⩊ Hoffman(1)^{⩊mi_d} を計算
                #mi_t = shuffle_product(HoffmanWordtoHoffman(mi_r),shpw(Hoffman(1),mi_d))
                mi_t = shuffle_product(HoffmanWordtoHoffman(mi_r),sh_y_pow(mi_d))

                # 係数を合わせる
                diffc = mi.coeff//mi_t.terms[mi.word]
                
                # zからひとつづつ引く
                z -= mi_t * diffc

                # r(正規化多項式)に足しこむ分
                zadd.terms[mi_r] = diffc
            end
            r.terms[mi_d] = zadd
        end
        return r
    else
        if w[1] != 2
            return Poly(HoffmanWordtoHoffman(w))
        end
        z = HoffmanWordtoHoffman(w)
        r = Poly{Hoffman}()
        
        while !iszero(z)
            m = leftmostys(z)
            zadd = Hoffman()
            mi_d = 0
            for mi in m
                mi_r, mi_d = leftysplit(mi.word)

                #mi_t = shuffle_product(HoffmanWordtoHoffman(mi_r),shpw(Hoffman(1),mi_d))
                mi_t = shuffle_product(HoffmanWordtoHoffman(mi_r),sh_y_pow(mi_d))
                
                diffc = mi.coeff//mi_t.terms[mi.word]

                z -= mi_t * diffc

                zadd.terms[mi_r] = diffc
            end
            r.terms[mi_d] = zadd
        end
        return r
    end
end
function shuffle_regularization_polynomial(h::Hoffman)::Poly{Hoffman}
    rp = Poly{Hoffman}()
    for (w,c) in h.terms
        rp += c * shuffle_regularization_polynomial(w)
    end
    return rp
end

# ρ

# 整数の配列とその配列の数の総和
struct CombSum
    a::Vector{Int}
    s::Int
end

"""
    sum_combinations(n,k)::Vector{CombSum}
和がn以下になるような2以上の整数k個の組をすべて返す
"""
function sum_combinations(n::UInt64, k::UInt64)::Vector{CombSum}

    if k == 0
        return CombSum[CombSum(Int[],0)]
    end

    results = Vector{CombSum}()

    a = fill(2, k)
    S = 2k  # 現在の総和

    while true
        push!(results, CombSum(copy(a), S))

        advanced = false
        suffix_sum = 0

        for i in k:-1:1
            suffix_sum += a[i]
            ai = a[i] + 1
            rest = k - i

            new_sum = (S - suffix_sum) + ai * (rest + 1)

            if new_sum <= n
                S = new_sum
                a[i] = ai
                for j in i+1:k
                    a[j] = ai
                end
                advanced = true
                break
            end
        end

        advanced || break
    end

    return results
end

"""
`Z^⧢ = ρ(Z^∗)` を満たす写像の`ρ(T^m)`を返す

Dict( (Tの次数,ゼータ関数の整数値の積) => 係数 )
という形で返す。

```jldoctest
julia> rho_t(UInt(2))
Dict{Tuple{Int64, Vector{Int64}}, Rational{BigInt}} with 2 entries:
  (2, [])  => 1
  (0, [2]) => 1
```

これは`ρ(T^2) = T^2 + ζ(2)`を表す

```jldoctest
julia> rho_t(UInt(4))
Dict{Tuple{Int64, Vector{Int64}}, Rational{BigInt}} with 5 entries:
  (4, [])     => 1
  (2, [2])    => 6
  (1, [3])    => -8
  (0, [2, 2]) => 3
  (0, [4])    => 6
```

これは`ρ(T^4) = T^4 + 6ζ(2)T^2 - 8ζ(3)T + 6ζ(4) + 3ζ(2)^2`を表す

多重ゼータ値は現れない

"""
function rho_t(m::UInt64)::Dict{Tuple{Int, Vector{Int}}, Rational{BigInt}}
    # key = (j, [n1,n2,...])  ↦  coefficient
    result = Dict{Tuple{Int, Vector{Int}}, Rational{BigInt}}()
    mfact = factorial(BigInt(m))

    for k in UInt64(0):m>>1
        m_k = mfact//factorial(BigInt(k))
        for cs in sum_combinations(m, k)
            #s = cs.s
            j = m - cs.s
            j < 0 && continue

            zetas = cs.a  # すでに nondecreasing
            coeff = m_k // factorial(BigInt(j))

            for n in zetas
                coeff *= (-1)^n // n
            end

            key = (j, zetas)
            result[key] = get(result, key, 0//1) + coeff
        end
    end

    return result
end
"""
`Z^⧢ = ρ(Z^∗)` を満たす写像
"""
function rho(ri::Poly{Index})::Poly{Hoffman}
    rh = Poly{Hoffman}()
    for (deg,coeff) in ri.terms
        # deg::Int, coeff::Index
        # c*T^d -> c*(rT)
        coeff_h = Hoffman(coeff)
        rT = rho_t(UInt64(deg))
        for (rhodegzeta,ratcoeff) in rT
            rhodeg = rhodegzeta[1]
            rhocoeff = coeff_h
            for rc in rhodegzeta[2]
                rhocoeff = shuffle_product(rhocoeff,IndexWordtoHoffman(Word(rc)))
            end
            rh.terms[rhodeg] = get(rh.terms, rhodeg, zero(Hoffman)) + rhocoeff
        end
    end
    return rh
end

# reg_st reg_sh
function reg_st(w::Word)::Index
    rp = stuffle_regularization_polynomial(w)
    return rp.terms[0]
end
function reg_st(idx::Index)::Index
    ri = Index()
    for (w,c) in idx.terms
        ri += c * reg_st(w)
    end
    return ri
end
function reg_sh(w::Word)::Hoffman
    rp = shuffle_regularization_polynomial(w)
    return rp.terms[0]
end
function reg_sh(h::Hoffman)::Hoffman
    rh = Hoffman()
    for (w,c) in h.terms
        rh += c * reg_sh(w)
    end
    return rh
end
