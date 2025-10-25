#[ monomials.jl ]#

# This file defines functions related to monomials

#=
export monomial_sh, monomial_st, monomial_st_star, monomial_sh_double, monomial_st_double, monomial_st_star_double,
       shuffle_product, stuffle_product, star_stuffle_product, 
       shuffle_product_double, stuffle_product_double, star_stuffle_product_double,
       shuffle_pow, stuffle_pow, star_stuffle_pow, 
       shpw, stpw, starstpw,
       Hoffman_hom, Hoffman_antihom, monomial_sw_w, starword_to_word,
       monomial_dual_h, monomial_dual_i, dual,
       monomial_hof_dual_h, monomial_hof_dual_i, Hoffman_dual,
       monomial_Landen, Landen_dual
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


###################################################################################################
############## products ###########################################################################

function shuffle_product(a::Hoffman, b::Hoffman)::Hoffman
    s1 = Hoffman()
    for (ma,ca) in a.terms
        s2 = Hoffman()
        for (mb,cb) in b.terms
            add!(s2,Hoffman(monomial_sh(ma.tovec,mb.tovec)),cb)
        end
        add!(s1,s2,ca)
    end
    return s1
end
function stuffle_product(a::Index, b::Index)::Index
    s1 = Index()
    for (ma,ca) in a.terms
        s2 = Index()
        for (mb,cb) in b.terms
            add!(s2,Index(monomial_st(ma.tovec,mb.tovec)),cb)
        end
        add!(s1,s2,ca)
    end
    return s1
end
function star_stuffle_product(a::Index, b::Index)::Index
    s1 = Index()
    for (ma,ca) in a.terms
        s2 = Index()
        for (mb,cb) in b.terms
            add!(s2,Index(monomial_st_star(ma.tovec,mb.tovec)),cb)
        end
        add!(s1,s2,ca)
    end
    return s1
end
function shuffle_product_double(a::Hoffman)::Hoffman
    result = Hoffman()
    pairs = collect(a.terms)
    n = lastindex(pairs)

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

    return result
end
function stuffle_product_double(a::Index)::Index
    result = Index()
    pairs = collect(a.terms)
    n = lastindex(pairs)

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

    return result
end
function star_stuffle_product_double(a::Index)::Index
    result = Index()
    pairs = collect(a.terms)
    n = lastindex(pairs)

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
        return shuffle_product(sp(a,n-1),a)
    end
end
function stpw(a::Index,n::Int)::Index
    r = one(Index)
    if n==0
        return r
    else
        return stuffle_product(sp(a,n-1),a)
    end
end
function starstpw(a::Index,n::Int)::Index
    r = one(Index)
    if n==0
        return r
    else
        return star_stuffle_product(sp(a,n-1),a)
    end
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

function starword_to_word(w::Hoffman)::Hoffman
    s = Hoffman()
    for (mw,cw) in w.terms
        add!(s,monomial_sw_w(mw),cw)
    end
    return s
end
starword_to_word(i::Index)::Index = Index(starword_to_word(i.toHoffman))


###################################################################################################
############## dual index #########################################################################


# 双対インデックス
function monomial_dual_h(w::Word)::Word
    r = Word((w[i]==1 ? 2 : 1) for i in lastindex(w):-1:1)
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
    return s
end

# Hoffman双対インデックス
function monomial_hof_dual_h(w::Word)::Word
    lw = lastindex(w)
    if lw == 0
        return w
    end
    return Word(2) * Tuple{Vararg{Int}}(3-i for i in w[2:lw])
end
function Hoffman_dual(w::Hoffman)::Hoffman
    s = Hoffman()
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
    c2 = -HoffmanWordtoHoffman(y)
    for i in 1:lastindex(w)
        if w[i] == 1
            s = s * c1
        else
            s = s * c2
        end
    end
    return s
end
function Landen_dual(w::Hoffman)::Hoffman
    s = Hoffman()
    for (mw,cw) in w.terms
        add!(s,monomial_Landen(mw),cw)
    end
    return s
end
Landen_dual(idx::Word)::Word = Index(monomial_Landen(IndexWordtoHoffmanWord(idx)))
Landen_dual(idx::Index)::Index = Index(Landen_dual(idx.toHoffman))

