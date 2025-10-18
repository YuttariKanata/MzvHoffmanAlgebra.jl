#[ monomials.jl ]#

# This file defines functions related to monomials

#=
export monomial_sh, monomial_st, monomial_sh_double, monomial_st_double,
       shuffle_product, stuffle_product, shuffle_product_double, stuffle_product_double,
       shuffle_pow, stuffle_pow, shpw,
       Hoffman_hom, Hoffman_antihom, monomial_sw_w, starword_to_word,
       monomial_dual, dual
=#

"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""

###################################################################################################
############## monomial functions #################################################################

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
function stuffle_product(a::Index, b::Index)
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


###################################################################################################
############## homomorphic ########################################################################

function Hoffman_hom(w::Word, image::Vector{Hoffman})::Hoffman
    s = one(Hoffman)
    for idx in w
        s *= image[idx]
    end
    return s
end
function Hoffman_antihom(w::Word, image::Vector{Hoffman})::Hoffman
    s = one(Hoffman)
    for idx in w
        s = image[idx] * s
    end
    return s
end

#star-valueのindexに対応するwordから通常のindexに
function monomial_sw_w(w::Word)::Hoffman
    return y * Hoffman_hom(w[2:end],[HoffmanWordtoHoffman(x),x + y])
end
function starword_to_word(w::Hoffman)::Hoffman
    s = Hoffman()
    for (mw,cw) in w.terms
        add!(s,monomial_sw_w(mw),cw)
    end
    return s
end

# 双対インデックス
function monomial_dual(w::Word)::Hoffman
    return Hoffman_antihom(w,[HoffmanWordtoHoffman(y),HoffmanWordtoHoffman(x)])
end
function dual(w::Hoffman)::Hoffman
    s = Hoffman()
    for (mw,cw) in w.terms
        add!(s,monomial_dual(mw),cw)
    end
    return s
end