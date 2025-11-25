#[ arithmetic.jl ]#

# This file defines arithmetic functions

import Base: +, -, *, ^, //

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

# isdefined table
#=

  Addition and Subtraction
-------------------------------------------------------------------------------------------------------------------
|  Add/Sub       |  NN  |  Word  |  Hoffman  |  MonoIndex  |  Index  |  Poly NN   |  Poly Hoffman  |  Poly Index  |
|----------------+------+--------+-----------+-------------+---------+------------+----------------+--------------|
|  NN            |  NN  |  Hof   |  Hof      |  Idx        |  Idx    |  Poly NN   |  Poly Hof      |  Poly Idx    |
|----------------+------+--------+-----------+-------------+---------+------------+----------------+--------------|
|  Word          |      |  Hof   |  Hof      |  X          |  X      |  Poly Hof  |  Poly Hof      |  X           |
|----------------+------+--------+-----------+-------------+---------+------------+----------------+--------------|
|  Hoffman       |      |        |  Hof      |  X          |  X      |  Poly Hof  |  Poly Hof      |  X           |
|----------------+------+--------+-----------+-------------+---------+------------+----------------+--------------|
|  MonoIndex     |      |        |           |  Idx        |  Idx    |  Poly Idx  |  X             |  Poly Idx    |
|----------------+------+--------+-----------+-------------+---------+------------+----------------+--------------|
|  Index         |      |        |           |             |  Idx    |  Poly Idx  |  X             |  Poly Idx    |
|----------------+------+--------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly NN       |      |        |           |             |         |  Poly NN   |  Poly Hof      |  Poly Idx    |
|----------------+------+--------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly Hoffman  |      |        |           |             |         |            |  Poly Hof      |  X           |
|----------------+------+--------+-----------+-------------+---------+------------+----------------+--------------|
|  Poly Index    |      |        |           |             |         |            |                |  Poly Idx    |
-------------------------------------------------------------------------------------------------------------------

=#



############################## ADDITIVE INVERSE ##############################

# Word
function -(a::Word)::Hoffman
    r = Hoffman()
    r.terms[a] = Rational(BigInt(-1))
    return r
end

# Hoffman
function -(a::Hoffman)::Hoffman
    result = Hoffman()
    for (w, c) in a.terms
        result.terms[w] = -c
    end
    return result
end

# MonoIndex
function -(a::MonoIndex)::MonoIndex
    result = MonoIndex(a.word,-a.coeff)
    return result
end

# Index
function -(a::Index)::Index
    result = Index()
    for (w, c) in a.terms
        result.terms[w] = -c
    end
    return result
end

# Poly
function -(a::Poly{A})::Poly{A} where A
    r = copy(a)
    for (d, h) in a.terms
        r.terms[d] = -h
    end
    return r
end



############################## ADD and SUBTRACT ##############################

########## Word ##########

# Word NN
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



########## Hoffman ##########

# Hoffman NN
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

# Hoffman Word
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



########## MonoIndex ##########

# MonoIndex NN
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

# MonoIndex MonoIndex
function +(a::MonoIndex, b::MonoIndex)::Index
    r = Index(a)
    if haskey(r.terms,b.word)
        if r.terms[b.word] == -b.coeff
            delete!(r.terms,b.word)
        else
            r.terms[b.word] += b.coeff
        end
    else
        r.terms[b.word] = b.coeff
    end
    return r
end
function -(a::MonoIndex, b::MonoIndex)::Index
    r = Index(a)
    if haskey(r.terms,b.word)
        if r.terms[b.word] == b.coeff
            delete!(r.terms,b.word)
        else
            r.terms[b.word] -= b.coeff
        end
    else
        r.terms[b.word] = -b.coeff
    end
    return r
end



########## Index ##########

# Index NN
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

# Index MonoIndex
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



########## Poly ##########

# auxiliary function
@inline function add!(r::Poly{A},b::ZetaExpr)::Poly{A} where A<:Union{NN,Hoffman,Index}
    if haskey(r.terms,0)
        if r.terms[0] == -A(b)
            delete!(r.terms,0)
        else
            r.terms[0] += b
        end
    else
        r.terms[0] = A(b)
    end
    return r
end
@inline function subtract!(r::Poly{A},b::ZetaExpr)::Poly{A} where A<:Union{NN,Hoffman,Index}
    if haskey(r.terms,0)
        if r.terms[0] == A(b)
            delete!(r.terms,0)
        else
            r.terms[0] -= b
        end
    else
        r.terms[0] = -A(b)
    end
    return r
end

########## Poly NN ##########

# Poly NN , NN
function +(a::Poly{Rational{BigInt}}, b::NN)::Poly{Rational{BigInt}}
    r = copy(a)
    return add!(r,b)
end
function -(a::Poly{Rational{BigInt}}, b::NN)::Poly{Rational{BigInt}}
    r = copy(a)
    return subtract!(r,b)
end

# Poly NN , (Word, Hoffman)
function +(a::Poly{Rational{BigInt}}, b::Union{Word,Hoffman})::Poly{Hoffman}
    r = Hoffman(a)
    return add!(r,b)
end
function -(a::Poly{Rational{BigInt}}, b::Union{Word,Hoffman})::Poly{Hoffman}
    r = Hoffman(a)
    return subtract!(r,b)
end

# Poly NN , (MonoIndex, Index)
function +(a::Poly{Rational{BigInt}}, b::Union{MonoIndex,Index})::Poly{Index}
    r = Index(a)
    return add!(r,b)
end
function -(a::Poly{Rational{BigInt}}, b::Union{MonoIndex,Index})::Poly{Index}
    r = Index(a)
    return subtract!(r,b)
end

# Poly NN , Poly NN
function +(a::Poly{Rational{BigInt}}, b::Poly{Rational{BigInt}})::Poly{Rational{BigInt}}
    r = copy(a)
    for (d,c) in b.terms
        if haskey(r.terms,d)
            if r.terms[d] == -c
                delete!(r.terms,d)
            else
                r.terms[d] += c
            end
            r.terms[d] = c
        end
    end
    return r
end
function -(a::Poly{Rational{BigInt}}, b::Poly{Rational{BigInt}})::Poly{Rational{BigInt}}
    r = copy(a)
    for (d,c) in b.terms
        if haskey(r.terms,d)
            if r.terms[d] == c
                delete!(r.terms,d)
            else
                r.terms[d] -= c
            end
            r.terms[d] = -c
        end
    end
    return r
end

+(a::Union{NN,Poly{Rational{BigInt}}},b::Poly{Rational{BigInt}})::Poly{Rational{BigInt}} = +(b,a)
-(a::Union{NN,Poly{Rational{BigInt}}},b::Poly{Rational{BigInt}})::Poly{Rational{BigInt}} = +(-b,a)
+(a::Union{Word,Hoffman},b::Poly{Rational{BigInt}})::Poly{Hoffman} = +(b,a)
-(a::Union{Word,Hoffman},b::Poly{Rational{BigInt}})::Poly{Hoffman} = +(-b,a)
+(a::Union{MonoIndex,Index},b::Poly{Rational{BigInt}})::Poly{Index} = +(b,a)
-(a::Union{MonoIndex,Index},b::Poly{Rational{BigInt}})::Poly{Index} = +(-b,a)



########## Poly Hoffman ##########

# Poly Hof , (NN, Word, Hoffman)
function +(a::Poly{Hoffman}, b::Union{NN,Word,Hoffman})::Poly{Hoffman}
    r = copy(a)
    return add!(r,b)
end
function -(a::Poly{Hoffman}, b::Union{NN,Word,Hoffman})::Poly{Hoffman}
    r = copy(a)
    return subtract!(r,b)
end

# Poly Hof , (Poly NN, Poly Hof)
function +(a::Poly{Hoffman}, b::Union{Poly{Rational{BigInt}},Poly{Hoffman}})::Poly{Hoffman}
    r = copy(a)
    for (d,c) in b.terms
        if haskey(r.terms,d)
            if r.terms[d] == -Hoffman(c)
                delete!(r.terms,d)
            else
                r.terms[d] += c
            end
            r.terms[d] = Hoffman(c)
        end
    end
    return r
end
function -(a::Poly{Hoffman}, b::Union{Poly{Rational{BigInt}},Poly{Hoffman}})::Poly{Hoffman}
    r = copy(a)
    for (d,c) in b.terms
        if haskey(r.terms,d)
            if r.terms[d] == Hoffman(c)
                delete!(r.terms,d)
            else
                r.terms[d] -= c
            end
            r.terms[d] = -Hoffman(c)
        end
    end
    return r
end

+(a::Union{NN,Word,Hoffman}, b::Poly{Hoffman})::Poly{Hoffman} = +(b,a)
-(a::Union{NN,Word,Hoffman}, b::Poly{Hoffman})::Poly{Hoffman} = +(-b,a)
+(a::Union{Poly{Rational{BigInt}},Poly{Hoffman}}, b::Poly{Hoffman})::Poly{Hoffman} = +(b,a)
-(a::Union{Poly{Rational{BigInt}},Poly{Hoffman}}, b::Poly{Hoffman})::Poly{Hoffman} = +(-b,a)



########## Poly Index ##########

# Poly{Index} , (NN,MonoIndex,Index)
function +(a::Poly{Index},b::Union{NN,MonoIndex,Index})::Poly{Index}
    r = copy(a)
    return add!(r,b)
end
function -(a::Poly{Index},b::Union{NN,MonoIndex,Index})::Poly{Index}
    r = copy(a)
    return subtract!(r,b)
end

# Poly{Index} , (Poly{Rational{BigInt}}, Poly{Index})
function +(a::Poly{Index},b::Union{Poly{Rational{BigInt}},Poly{Index}})::Poly{Index}
    r = copy(a)
    for (d,c) in b.terms
        if haskey(r.terms,d)
            if r.terms[d] == -Index(c)
                delete!(r.terms,d)
            else
                r.terms[d] += c
            end
            r.terms[d] = Index(c)
        end
    end
    return r
end
function -(a::Poly{Index},b::Union{Poly{Rational{BigInt}},Poly{Index}})::Poly{Index}
    r = copy(a)
    for (d,c) in b.terms
        if haskey(r.terms,d)
            if r.terms[d] == Index(c)
                delete!(r.terms,d)
            else
                r.terms[d] -= c
            end
            r.terms[d] = -Index(c)
        end
    end
    return r
end

+(a::Union{NN,MonoIndex,Index}, b::Poly{Index})::Poly{Index} = +(b,a)
-(a::Union{NN,MonoIndex,Index}, b::Poly{Index})::Poly{Index} = +(-b,a)
+(a::Union{Poly{Rational{BigInt}},Poly{Index}}, b::Poly{Index})::Poly{Index} = +(b,a)
-(a::Union{Poly{Rational{BigInt}},Poly{Index}}, b::Poly{Index})::Poly{Index} = +(-b,a)





############################## MULTIPLICATION ##############################
########## Word ##########
# Word NN
function *(a::NN, b::Word)::Hoffman
    r = Hoffman()
    if a != 0
        r.terms[b] = a
    end
    return r
end
*(a::Word, b::NN)::Hoffman = *(b,a)

# Word Word
@inline function *(a::Word, b::Word)::Word
    return Word(a... , b...)
end



########## Hoffman ##########
# Hoffman NN
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

# Hoffman Word
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



########## MonoIndex ##########
# MonoIndex NN
function *(a::NN, b::MonoIndex)::Index
    r = Index()
    if a != 0
        r.terms[b.word] = a * b.coeff
    end
    return r
end
*(a::MonoIndex, b::NN)::Index = *(b,a)

# MonoIndex MonoIndex
function *(a::MonoIndex, b::MonoIndex)::MonoIndex
    return MonoIndex(a.word*b.word,a.coeff*b.coeff)
end



########## Index ##########
# Index NN
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

# Index MonoIndex
function *(a::Index, b::MonoIndex)::Index
    r = Index()
    for (w,c) in a.terms
        r.terms[w*b.word] = c*b.coeff
    end
    return r
end
function *(a::MonoIndex, b::Index)::Index
    r = Index()
    for (w,c) in b.terms
        r.terms[a.word*w] = c*a.coeff
    end
    return r
end

# Index Index
function *(a::Index, b::Index)::Index
    result = Index()
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



########## Poly NN ##########
# PolyNN , NN
function *(a::NN, b::Poly{Rational{BigInt}})::Poly{Rational{BigInt}}
    r = Poly{Rational{BigInt}}()
    if a != 0
        for (d, h) in b.terms
            r.terms[d] = a * h
        end
    end
    return r
end
*(a::Poly{Rational{BigInt}}, b::NN)::Poly{Rational{BigInt}} = *(b,a)

# Poly NN , (Word, Hoffman)
function *(a::Poly{Rational{BigInt}}, b::Union{Word,Hoffman})::Poly{Hoffman}
    r = Poly{Hoffman}()
    for (d,c) in a.terms
        r.terms[d] = c*b
    end
    return r
end
*(a::Union{Word,Hoffman},b::Poly{Rational{BigInt}})::Poly{Hoffman} = *(b,a)

# Poly NN , (MonoIndex, Index)
function *(a::Poly{Rational{BigInt}}, b::Union{MonoIndex,Index})::Poly{Index}
    r = Poly{Index}()
    for (d,c) in a.terms
        r.terms[d] = c*b
    end
    return r
end
*(a::Union{MonoIndex,Index},b::Poly{Rational{BigInt}})::Poly{Index} = *(b,a)

# Poly NN , Poly NN
function normal_multiply_Poly(a::Poly{Rational{BigInt}}, b::Poly{Rational{BigInt}})::Poly{Rational{BigInt}}
    r = Poly{Rational{BigInt}}()
    for (da,ha) in a.terms
        for (db,hb) in b.terms
            deg = da + db    # 次数の足し算
            prod = ha * hb   # 係数の掛け算
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
function *(a::Poly{Rational{BigInt}}, b::Poly{Rational{BigInt}})::Poly{Rational{BigInt}}
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
    return normal_multiply_Poly(a,b)
end



########## Poly Hoffman ##########
# Poly Hoffman, NN
function *(a::NN, b::Poly{Hoffman})::Poly{Hoffman}
    r = Poly{Hoffman}()
    if a != 0
        for (d, h) in b.terms
            r.terms[d] = a * h
        end
    end
    return r
end
*(a::Poly{Hoffman}, b::NN)::Poly{Hoffman} = *(b,a)

# Poly Hoffman , (Word, Hoffman)
function *(a::Poly{Hoffman}, b::Union{Word,Hoffman})::Poly{Hoffman}
    r = Poly{Hoffman}()
    for (d,c) in a.terms
        r.terms[d] = c*b
    end
    return r
end
function *(a::Union{Word,Hoffman}, b::Poly{Hoffman})::Poly{Hoffman}
    r = Poly{Hoffman}()
    for (d,c) in b.terms
        r.terms[d] = a*c
    end
    return r
end

# Poly Hoffman , Poly NN
function normal_multiply_Poly(a::Poly{Hoffman}, b::Poly{Rational{BigInt}})::Poly{Hoffman}
    r = Poly{Hoffman}()
    for (da,ha) in a.terms
        for (db,hb) in b.terms
            deg = da + db    # 次数の足し算
            prod = ha * hb   # 係数の掛け算
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
normal_multiply_Poly(a::Poly{Rational{BigInt}}, b::Poly{Hoffman})::Poly{Hoffman} = normal_multiply_Poly(b,a)
function *(a::Poly{Hoffman}, b::Poly{Rational{BigInt}})::Poly{Hoffman}
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
    return normal_multiply_Poly(a,b)
end
*(a::Poly{Rational{BigInt}}, b::Poly{Hoffman})::Poly{Hoffman} = *(b,a)

# Poly Hoffman , Poly Hoffman
function normal_multiply_Poly(a::Poly{Hoffman}, b::Poly{Hoffman})::Poly{Hoffman}
    r = Poly{Hoffman}()
    for (da,ha) in a.terms
        for (db,hb) in b.terms
            deg = da + db    # 次数の足し算
            prod = ha * hb   # 係数の掛け算
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
function *(a::Poly{Hoffman}, b::Poly{Hoffman})::Poly{Hoffman}
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
    return normal_multiply_Poly(a,b)
end



########## Poly Index ##########
# Poly Index, NN
function *(a::NN, b::Poly{Index})::Poly{Index}
    r = Poly{Index}()
    if a != 0
        for (d, i) in b.terms
            r.terms[d] = a * i
        end
    end
    return r
end
*(a::Poly{Index}, b::NN)::Poly{Index} = *(b,a)

# Poly Index , (MonoIndex, Index)
function *(a::Poly{Index}, b::Union{MonoIndex,Index})::Poly{Index}
    r = Poly{Index}()
    for (d,c) in a.terms
        r.terms[d] = c*b
    end
    return r
end
function *(a::Union{MonoIndex,Index}, b::Poly{Index})::Poly{Index}
    r = Poly{Index}()
    for (d,c) in b.terms
        r.terms[d] = a*c
    end
    return r
end

# Poly Index , Poly NN
function normal_multiply_Poly(a::Poly{Index}, b::Poly{Rational{BigInt}})::Poly{Index}
    r = Poly{Index}()
    for (da,ha) in a.terms
        for (db,hb) in b.terms
            deg = da + db    # 次数の足し算
            prod = ha * hb   # 係数の掛け算
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
normal_multiply_Poly(a::Poly{Rational{BigInt}}, b::Poly{Index})::Poly{Index} = normal_multiply_Poly(b,a)
function *(a::Poly{Index}, b::Poly{Rational{BigInt}})::Poly{Index}
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
    return normal_multiply_Poly(a,b)
end
*(a::Poly{Rational{BigInt}}, b::Poly{Index})::Poly{Index} = *(b,a)

# Poly Index , Poly Index
function normal_multiply_Poly(a::Poly{Index}, b::Poly{Index})::Poly{Index}
    r = Poly{Index}()
    for (da,ha) in a.terms
        for (db,hb) in b.terms
            deg = da + db    # 次数の足し算
            prod = ha * hb   # 係数の掛け算
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
function *(a::Poly{Index}, b::Poly{Index})::Poly{Index}
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
    return normal_multiply_Poly(a,b)
end



############################## POWER ##############################
# Word
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

# Hoffman
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

# MonoIndex
function ^(a::MonoIndex, n::Integer)::MonoIndex
    return MonoIndex(a.word^n,a.coeff^n)
end

# Index
function ^(a::Index, n::Integer)::Index
    if n < 0
        throw(ArgumentError("Index-type powers for negative exponents are not defined"))
    elseif n == 0
        return one(Index)
    elseif n == 1
        return copy(a)
    elseif isone(a)
        return one(Index)
    end

    result = one(Index)
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

# Poly
function ^(a::Poly{A}, n::Integer)::Poly{A} where A
    if n < 0
        throw(DomainError(n, "negative power is not supported for Poly"))
    elseif n == 0
        return one(Poly{A})
    elseif n == 1
        return copy(a)
    end

    # 最適化: a が単項 T^k (係数が exactly one(Hoffman)) の場合
    if length(a.terms) == 1
        (deg, coeff) = first(a.terms)
        r = Poly{A}()
        if deg == 0
            r.terms[0] = coeff ^ n
        elseif isone(coeff)
            r.terms[deg * Int(n)] = one(A)
        else
            r.terms[deg * Int(n)] = coeff ^ n
        end
        return r
    end    

    # 二乗累乗（binary exponentiation）
    base = copy(a)
    result = one(Poly{A})
    nn = Int(n)  # n は非負なので安全に Int に落とす

    while nn > 0
        if (nn & 1) == 1
            result = result * base   # 既に定義済みの *(Poly, Poly) を使用
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
function shift_degree(r::Poly{A},n::Int)::Poly{A} where A
    out = Poly{A}()
    for (d, h) in r.terms
        out.terms[d+n] = copy(h)
    end
    return out
end

# auxiliary function
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



############################# DIVISION for NN ##############################
# NN Word
function //(b::Word, a::NN)::Hoffman
    r = Hoffman()
    r.terms[b] = 1//a
    return r
end

# NN Hoffman
function //(b::Hoffman, a::NN)::Hoffman
    r = Hoffman()
    for (w, c) in b.terms
        r.terms[w] = c//a
    end
    return r
end

# NN MonoIndex
function //(b::MonoIndex, a::NN)::Index
    r = Index()
    r.terms[b.word] = b.coeff//a
    return r
end

# NN Index
function //(b::Index, a::NN)::Index
    r = Index()
    for (w, c) in b.terms
        r.terms[w] = c//a
    end
    return r
end

# NN Poly
function //(b::Poly{A}, a::NN)::Poly{A} where A
    r = Poly{A}()
    for (d, h) in b.terms
        r.terms[d] = h // a
    end 
    return r
end