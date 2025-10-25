#[ arithmetic.jl ]#

# This file defines arithmetic functions

import Base: +, -, *, ^

#=
export shift_degree, add!
=#

"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## Arithmetic Functions ###############################################################


############################## ADD and SUBTRACT ##############################
########## NN ##########
# NN Word
function +(a::Word, b::NN)::Hoffman
    w = Hoffman()
    if a != one(Word)      # one(Word)
        w.terms[a] = Rational(BigInt(1))
        w.terms[Word()] = b
    else
        if b != Clong(-1)
            w.terms[Word()] = b + 1
        end
    end
    return w
end
+(a::NN, b::Word)::Hoffman = +(b,a)
-(a::Word, b::NN)::Hoffman = +(a,-b)
function -(a::NN, b::Word)::Hoffman
    w = Hoffman()
    if b != one(Word)      # one(Word)
        w.terms[b] = Rational(BigInt(-1))
        w.terms[Word()] = a
    else
        if a != Culong(1)
            w.terms[Word()] = a - 1
        end
    end
    return w
end

# NN Hoffman
function +(a::Hoffman, b::NN)::Hoffman
    r = copy(a)
    if haskey(r.terms,Word())
        if r.terms[Word()] == -b
            delete!(r.terms,Word())
        else
            r.terms[Word()] += b
        end
    else
        r.terms[Word()] = b
    end
    return r
end
+(a::NN, b::Hoffman)::Hoffman = +(b,a)
-(a::Hoffman, b::NN)::Hoffman = +(a,-b)
-(a::NN, b::Hoffman)::Hoffman = +(-b,a)

# NN Index
function +(a::Index, b::NN)::Index
    result = copy(a)
    if haskey(result.terms,Word())
        if result.terms[Word()] == -b
            delete!(result.terms,Word())
        else
            result.terms[Word()] += b
        end
    else
        result.terms[Word()] = b
    end
    return result
end
+(a::NN, b::Index)::Index = +(b,a)
-(a::Index, b::NN)::Index = +(a,-b)
-(a::NN, b::Index)::Index = +(-b,a)

# NN MonoIndex
function +(a::MonoIndex, b::NN)::Index
    r = Index()
    if a.word != Word()
        r.terms[a.word] = Rational(BigInt(1))
        r.terms[Word()] = b
    else
        if b != Clong(-1)
            r.terms[Word()] = b + 1
        end
    end
    return r
end
+(a::NN, b::MonoIndex)::Index = +(b,a)
-(a::MonoIndex, b::NN)::Index = +(a,-b)
-(a::NN, b::MonoIndex)::Index = +(-b,a)

## NN RegHoffman
function +(a::RegHoffman, b::NN)::RegHoffman
    r = copy(a)
    if haskey(r.terms,0)
        if r.terms[0] == Hoffman(-b)
            delete!(r.terms,0)
        else
            r.terms[0] += b
        end
    else
        r.terms[0] = b
    end
    return r
end
+(a::NN, b::RegHoffman)::RegHoffman = +(b,a)
-(a::RegHoffman, b::NN)::RegHoffman = +(a,-b)
-(a::NN, b::RegHoffman)::RegHoffman = +(-b,a)


########## Word ##########
# Word Word
function +(a::Word, b::Word)::Hoffman
    r = Hoffman()
    if a != b
        r.terms[a] = Rational(BigInt(1))
        r.terms[b] = Rational(BigInt(1))
        return r
    else
        r.terms[a] = Rational(BigInt(2))
        return r
    end
end
function -(a::Word, b::Word)::Hoffman
    r = Hoffman()
    if a != b
        r.terms[a] = Rational(BigInt(1))
        r.terms[b] = Rational(BigInt(-1))
    end # if a==b then return Hoffman()
    return r
end

# Word Hoffman
function +(a::Hoffman, b::Word)::Hoffman
    r = copy(a)
    if haskey(r.terms,b)
        if r.terms[b] == -1
            delete!(r.terms,b)
        else
            r.terms[b] += 1
        end
    else
        r.terms[b] = 1
    end
    return r
end
+(a::Word, b::Hoffman)::Hoffman = +(b,a)
function -(a::Hoffman, b::Word)::Hoffman
    r = copy(a)
    if haskey(r.terms,b)
        if r.terms[b] == Culong(1)
            delete!(r.terms,b)
        else
            r.terms[b] -= Rational(BigInt(1))
        end
    else
        r.terms[b] = Rational(BigInt(-1))
    end
    return r
end
-(a::Word, b::Hoffman)::Hoffman = +(-b,a)

# Word RegHoffman
function +(a::RegHoffman, b::Word)::RegHoffman
    r = copy(a)
    if haskey(r.terms,b)
        if r.terms[b] == Hoffman(-1)
            delete!(r.terms,b)
        else
            r.terms[b] += Hoffman(1)
        end
    else
        r.terms[b] = Hoffman(1)
    end
    return r
end
+(a::Word, b::RegHoffman)::RegHoffman = +(b,a)
function -(a::RegHoffman, b::Word)::RegHoffman
    r = copy(a)
    if haskey(r.terms,b)
        if r.terms[b] == Hoffman(1)
            delete!(r.terms,b)
        else
            r.terms[b] += Hoffman(-1)
        end
    else
        r.terms[b] = Hoffman(-1)
    end
    return r
end
-(a::Word, b::RegHoffman)::RegHoffman = +(-b,a)


########## Hoffman ##########
# Hoffman
function -(a::Hoffman)::Hoffman
    result = Hoffman()
    for (w, c) in a.terms
        result.terms[w] = -c
    end
    return result
end

# Hoffman Hoffman
function +(a::Hoffman, b::Hoffman)::Hoffman
    result = copy(a)
    # b の項をマージ
    for (w, c) in b.terms
        if haskey(result.terms, w)
            # 係数が 0 になったら削除
            if result.terms[w] == -c
                delete!(result.terms, w)
            else
                result.terms[w] += c
            end
        else
            result.terms[w] = c
        end
    end
    return result
end
-(a::Hoffman,b::Hoffman)::Hoffman = +(a,-b)

# Hoffman RegHoffman
function +(a::RegHoffman, b::Hoffman)::RegHoffman
    r = copy(a)
    if haskey(r.terms, 0)
        r.terms[0] += b
        if iszero(r.terms[0])
            delete!(r.terms, 0)
        end
    else
        r.terms[0] = copy(b)
    end
    return r
end
+(a::Hoffman, b::RegHoffman)::RegHoffman = +(b,a)
-(a::RegHoffman, b::Hoffman)::RegHoffman = +(a,-b)
-(a::Hoffman, b::RegHoffman)::RegHoffman = +(-b,a)


########## MonoIndex ##########
# MonoIndex
function -(a::MonoIndex)::MonoIndex
    result = MonoIndex(a.word,-a.coeff)
    return result
end

# MonoIndex Index
function +(a::Index, b::MonoIndex)::Index
    result = copy(a)
    if haskey(result.terms, b.word)
        # 係数が 0 になったら削除
        if result.terms[b.word] == -b.coeff
            delete!(result.terms, b.word)
        else
            result.terms[b.word] += b.coeff
        end
    else
        result.terms[b.word] = b.coeff
    end
    return result
end
+(a::MonoIndex,b::Index)::Index = +(b,a)
function -(a::Index, b::MonoIndex)::Index
    result = copy(a)
    if haskey(result.terms, b.word)
        # 係数が 0 になったら削除
        if result.terms[b.word] == b.coeff
            delete!(result.terms, b.word)
        else
            result.terms[b.word] -= b.coeff
        end
    else
        result.terms[b.word] = -b.coeff
    end
    return result
end
function -(a::MonoIndex, b::Index)::Index
    result = Index()
    for (w, c) in b.terms
        result.terms[w] = -c
    end
    if haskey(result.terms, a.word)
        # 係数が 0 になったら削除
        if result.terms[a.word] == -a.coeff
            delete!(result.terms, a.word)
        else
            result.terms[a.word] += a.coeff
        end
    else
        result.terms[a.word] = a.coeff
    end
    return result
end


########## Index ##########
# Index
function -(a::Index)::Index
    result = Index()
    for (w, c) in a.terms
        result.terms[w] = -c
    end
    return result
end

# Index Index
function +(a::Index, b::Index)::Index
    result = copy(a)
    # b の項をマージ
    for (w, c) in b.terms
        if haskey(result.terms, w)
            # 係数が 0 になったら削除
            if result.terms[w] == -c
                delete!(result.terms, w)
            else
                result.terms[w] += c
            end
        else
            result.terms[w] = c
        end
    end
    return result
end
function -(a::Index, b::Index)::Index
    result = copy(a)
    # b の項をマージ
    for (w, c) in b.terms
        if haskey(result.terms, w)
            # 係数が 0 になったら削除
            if result.terms[w] == c
                delete!(result.terms, w)
            else
                result.terms[w] -= c
            end
        else
            result.terms[w] = -c
        end
    end
    return result
end


########## RegHoffman ##########
# RegHoffman
function -(a::RegHoffman)::RegHoffman
    r = copy(a)
    for (d, h) in a.terms
        r.terms[d] = -h
    end
    return r
end

# RegHoffman RegHoffman
function +(a::RegHoffman,b::RegHoffman)::RegHoffman
    r = copy(a)
    for (deg,hb) in b.terms
        if haskey(r.terms, deg)
            s = r.terms[deg] + hb
            if iszero(s)
                delete!(r.terms, deg)
            else
                r.terms[deg] = s
            end
        else
            r.terms[deg] = copy(hb)
        end
    end
    return r
end
-(a::RegHoffman, b::RegHoffman)::RegHoffman = +(a,-b)




############################## MULTIPLICATION ##############################
########## NN ##########
# NN Word
function *(a::NN, b::Word)::Hoffman
    r = Hoffman()
    if a != 0
        r.terms[b] = a
    end
    return r
end
*(a::Word, b::NN)::Hoffman = *(b,a)

# NN Hoffman
function *(a::NN, b::Hoffman)::Hoffman
    r = Hoffman()
    if a != 0
        for (w, c) in b.terms
            r.terms[w] = c*a
        end
    end
    return r
end
*(a::Hoffman, b::NN)::Hoffman = *(b,a)

# NN MonoIndex
function *(a::NN, b::MonoIndex)::Index
    r = Index()
    if a != 0
        r.terms[b.word] = a * b.coeff
    end
    return r
end
*(a::MonoIndex, b::NN)::Index = *(b,a)

# NN Index
function *(a::NN, b::Index)::Index
    r = Index()
    if a != 0
        for (w, c) in b.terms
            r.terms[w] = c * a
        end
    end
    return r
end
*(a::Index, b::NN)::Index = *(b,a)

# NN RegHoffman
function *(a::NN, b::RegHoffman)::RegHoffman
    r = RegHoffman()
    if a != 0
        for (d, h) in b.terms
            r.terms[d] = a * h
        end
    end
    return r
end
*(a::RegHoffman, b::NN)::RegHoffman = *(b,a)


########## Word ##########
# Word Word
@inline function *(a::Word, b::Word)::Word
    return Word(a... , b...)
end

@inline function ^(t::Word, n::Integer)::Word
    n <= 0 && return Word()
    len = length(t)
    total = len * n
    v = Vector{eltype(t)}(undef, total)
    for i in 0:n-1
        copyto!(v, i*len + 1, t, 1, len)
    end
    return Word(v)
end


# Word Hoffman
function *(a::Word, b::Hoffman)::Hoffman
    result = Hoffman()
    for (wb,cb) in b.terms
        result.terms[a*wb] = cb
    end
    return result
end
function *(a::Hoffman,b::Word)::Hoffman
    result = Hoffman()
    for (wa,ca) in a.terms
        result.terms[wa*b] = ca
    end
    return result
end

# Word RegHoffman
function *(a::Word, b::RegHoffman)::RegHoffman
    r = RegHoffman()
    for (d, h) in b.terms
        r.terms[d] = a * h
    end
    return r
end
function *(a::RegHoffman, b::Word)::RegHoffman
    r = RegHoffman()
    for (d, h) in a.terms
        r.terms[d] = h * b
    end
    return r
end


########## Hoffman ##########
# Hoffman Hoffman
function *(a::Hoffman, b::Hoffman)::Hoffman
    result = Hoffman()
    for (wa, ca) in a.terms
        for (wb, cb) in b.terms
            wprod = wa * wb
            coeff = ca * cb
            if haskey(result.terms, wprod)
                result.terms[wprod] += coeff
            else
                result.terms[wprod] = coeff
            end
        end
    end
    filter!(p->!iszero(p.second),result.terms)
    return result
end
function ^(a::Hoffman, n::Integer)::Hoffman
    if n < 0
        throw(ArgumentError("Hoffman-type powers for negative exponents are not defined"))
    elseif n == 0
        return one(Hoffman)
    elseif n == 1
        return copy(a)
    elseif isone(a)
        return one(Hoffman)
    end

    result = one(Hoffman)
    base = copy(a)
    nn = n

    while nn > 0
        if (nn & 1) == 1
            result = result * base
        end
        nn >>= 1
        if nn > 0
            base = base * base
        end
    end

    return result
end


# Hoffman RegHoffman
function *(a::RegHoffman, b::Hoffman)::RegHoffman
    r = RegHoffman()
    for (deg, h) in a.terms
        r.terms[deg] = h * b
    end
    return r
end
function *(a::Hoffman, b::RegHoffman)::RegHoffman
    r = RegHoffman()
    for (deg, h) in b.terms
        r.terms[deg] = a * h
    end
    return r
end

########## RegHoffman ##########
# RegHoffman RegHoffman
function normal_multiply_RegHoffman(a::RegHoffman, b::RegHoffman)::RegHoffman
    r = RegHoffman()
    for (da,ha) in a.terms
        for (db,hb) in b.terms
            deg = da + db
            prod = ha * hb   # Hoffman の掛け算
            if haskey(r.terms, deg)
                r.terms[deg] += prod
            else
                r.terms[deg] = prod
            end
        end
    end
    filter!(p->!iszero(p.second),r.terms)
    return r
end
function *(a::RegHoffman, b::RegHoffman)::RegHoffman
    # a が T^n 型なら次数を足すだけ
    if length(a.terms) == 1 && isone(first(values(a.terms)))
        n = first(keys(a.terms))
        return shift_degree(b, n)
    end
    # b が T^n 型なら同様
    if length(b.terms) == 1 && isone(first(values(b.terms)))
        n = first(keys(b.terms))
        return shift_degree(a, n)
    end
    # 通常の掛け算
    return normal_multiply_RegHoffman(a,b)
end
function ^(a::RegHoffman, n::Integer)::RegHoffman
    if n < 0
        throw(DomainError(n, "negative power is not supported for RegHoffman"))
    elseif n == 0
        return one(RegHoffman)
    elseif n == 1
        return copy(a)
    end

    # 最適化: a が単項 T^k (係数が exactly one(Hoffman)) の場合
    if length(a.terms) == 1
        (deg, coeff) = first(a.terms)
        r = RegHoffman()
        if deg == 0
            r.terms[0] = coeff ^ n
        elseif isone(coeff)
            r.terms[deg * Int(n)] = one(Hoffman)
        else
            r.terms[deg * Int(n)] = coeff ^ n
        end
        return r
    end    

    # 二乗累乗（binary exponentiation）
    base = copy(a)
    result = one(RegHoffman)
    nn = Int(n)  # n は非負なので安全に Int に落とす

    while nn > 0
        if (nn & 1) == 1
            result = result * base   # 既に定義済みの *(RegHoffman, RegHoffman) を使用
        end
        nn >>= 1
        if nn > 0
            base = base * base
        end
    end

    return result
end

# degree shift
""" r , n -> r T^n """
function shift_degree(r::RegHoffman,n::Int)::RegHoffman
    out = RegHoffman()
    for (d, h) in r.terms
        out.terms[d+n] = copy(h)
    end
    return out
end





""" h += w*c """
function add!(h::Hoffman,w::Hoffman,c::Union{Rational,Integer} = Rational(BigInt(1)))

    for (ww,wc) in w.terms
        if haskey(h.terms,ww)
            if h.terms[ww] == -wc*c
                delete!(h.terms,ww)
            else
                h.terms[ww] += wc*c
            end
        else
            h.terms[ww] = wc*c
        end
    end
end
function add!(h::Index,w::Index,c::Union{Rational,Integer} = Rational(BigInt(1)))

    for (ww,wc) in w.terms
        if haskey(h.terms,ww)
            if h.terms[ww] == -wc*c
                delete!(h.terms,ww)
            else
                h.terms[ww] += wc*c
            end
        else
            h.terms[ww] = wc*c
        end
    end
end
