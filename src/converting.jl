#[ converting.jl ]#

# This file defines conversion functions between each type

#=
export HoffmanWordtoMonoIndex, IndexWordtoMonoIndex,
       word, IndexWordtoHoffmanWord, HoffmanWordtoIndexWord,
       HoffmanWordtoIndex, IndexWordtoIndex,
       HoffmanWordtoHoffman, IndexWordtoHoffman,
       x, y, T
=#

"""
###################################################################################################
                                        Operators
###################################################################################################
"""


###################################################################################################
############## Conversion Functions ###############################################################


function Operator(sym::AbstractString, n::Integer=1)::Operator
    op = Operator()
    push!(op.ops,str_to_op(sym,n))
    return op
end
function Operator(v::Vector{<:AbstractOp})::Operator
    opes = Operator()
    append_clean!(opes,v)
    return opes
end
function Operator(v::Vector{<:Tuple{AbstractString, Integer}})::Operator
    lv = lastindex(v)
    regv = Vector{AbstractOp}(undef,lv)
    for i in 1:lv
        regv[i] = str_to_op(v[i])
    end
    opes = Operator()
    append_clean!(opes,regv)
    return opes
end
Operator(op::Operator)::Operator = op
function Operator(op::AbstractOp)::Operator
    r = Operator()
    push!(r.ops,op)
    return r
end
function Operator(op::Type{<:AbstractOp})::Operator
    r = Operator()
    if op === OpDown
        push!(r.ops,OpUp(-1))
    else
        push!(r.ops,(op)(1))
    end
    return r
end
function Operator(op::Type{<:AbstractOp},n::Int,k::Int=1)::Operator
    r = Operator()
    if op === OpDown
        push!(r.ops,OpUp(-n))
    elseif op === OpDeriv
        push!(r.ops,OpDeriv(k,n))
    else
        push!(r.ops,(op)(n))
    end
    return r
end





"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## Conversion Functions ###############################################################

# [============== about MonoIndex ==============]
function MonoIndex(v::Union{Vector,Word})::MonoIndex
    return MonoIndex(Tuple{Vararg{ExprInt}}(v),Rational(BigInt(1)))
end
function MonoIndex(x::Index)::MonoIndex
    if !is_monomial(x)
        throw(DomainError(x,"x must be monomial"))
    end
    return MonoIndex(first(keys(x.terms)),first(values(x.terms)))
end
function MonoIndex(w::Hoffman)::MonoIndex
    if !is_monomial(w)
        throw(DomainError(w,"w must be monomial"))
    end
    wo = first(keys(w.terms))
    if wo[1] != 2 # yに対応するだけ
        throw(DomainError(w,"w does not start with y"))
    end
    return MonoIndex(idxprs(wo))
end

function HoffmanWordtoMonoIndex(w::Word)::MonoIndex
    return MonoIndex(idxprs(w))
end
function IndexWordtoMonoIndex(w::Word)::MonoIndex
    return MonoIndex(w,Rational(BigInt(1)))
end

MonoIndex(m::MonoIndex)::MonoIndex = m
MonoIndex(n::NN)::MonoIndex = MonoIndex((),n)
# TODO: MonoIndex for hrm shf mpl 


# [============== about Word ==============]
word(x::Int...)::Word = Word(x)
const x = word(1)
const y = word(2)

# Vector -> Word は Word(v) でいける
word(m::MonoIndex)::Word = isone(m.coeff) ? m.word : throw(DomainError(m,"w's coefficient is not 1"))
"""returned Index word"""
function word(i::Index)::Word
    if !is_monomial(i)
        throw(DomainError(i,"i must be monomial"))
    end
    if !isone(first(values(i.terms)))
        throw(DomainError(i,"i's coefficient is not 1"))
    end
    return Word(first(keys(i.terms)))
end
"""returned Hoffman word"""
function word(w::Hoffman)::Word
    if !is_monomial(w)
        throw(DomainError(i,"i must be monomial"))
    end
    if !isone(first(values(i.terms)))
        throw(DomainError(i,"i's coefficient is not 1"))
    end
    return Word(first(keys(i.terms)))
end

function IndexWordtoHoffmanWord(w::Word)::Word
    return idxdprs(collect(w))
end
function HoffmanWordtoIndexWord(w::Word)::Word
    return Word(idxprs(w))
end

word(w::Word)::Word = w
# TODO: Word for hrm shf mpl


# [============== about Index ==============]
function Index(v::Vector{Vector{Int}})::Index
    idx = Index()
    c = Rational(BigInt(1))
    for vi in v
        wvi = Word(vi)
        if haskey(idx.terms,wvi)
            idx.terms[wvi] += c
        else
            idx.terms[wvi] = c
        end
    end
    return idx
end
function Index(m::MonoIndex)::Index
    idx = Index()
    idx.terms[m.word] = m.coeff
    return idx
end
function HoffmanWordtoIndex(w::Word,coeff::NN = Rational(BigInt(1)))::Index
    idx = Index()
    idx.terms[HoffmanWordtoIndexWord(w)] = coeff
    return idx
end
function IndexWordtoIndex(w::Word,coeff::NN = Rational(BigInt(1)))::Index
    idx = Index()
    idx.terms[w] = coeff
    return idx
end
function Index(w::Hoffman)::Index
    idx = Index()
    for (k,v) in w.terms
        idx.terms[HoffmanWordtoIndexWord(k)] = v
    end
    return idx
end

Index(i::Index)::Index = i
function Index(n::NN)::Index
    idx = Index()
    idx.terms[()] = n
    return idx
end
# TODO: Index for hrm shf mpl


# [============== about Hoffman ==============]
function Hoffman(v::Vector{Vector{Int}})::Hoffman
    w = Hoffman()
    c = Rational(BigInt(1))
    for vi in v
        wvi = Word(vi)
        if haskey(w.terms,wvi)
            w.terms[wvi] += c
        else
            w.terms[wvi] = c
        end
    end
    return w
end
function Hoffman(m::MonoIndex)::Hoffman
    w = Hoffman()
    w.terms[IndexWordtoHoffmanWord(m.word)] = m.coeff
    return w
end
function HoffmanWordtoHoffman(wm::Word)::Hoffman
    w = Hoffman()
    w.terms[wm] = Rational(BigInt(1))
    return w
end
function IndexWordtoHoffman(i::Word)::Hoffman
    w = Hoffman()
    w.terms[IndexWordtoHoffmanWord(i)] = Rational(BigInt(1))
    return w
end
function Hoffman(i::Index)::Hoffman
    w = Hoffman()
    for (k,v) in i.terms
        w.terms[IndexWordtoHoffmanWord(k)] = v
    end
    return w
end

Hoffman(w::Hoffman)::Hoffman = w
function Hoffman(n::NN)::Hoffman
    w = Hoffman()
    w.terms[()] = n
    return w
end
# TODO: Hoffman for hrm shf mpl


# [============== about RegHoffman ==============]
const T = let
    t = RegHoffman()
    t.terms[1] = one(Hoffman)
    t
end
function RegHoffman(w::Hoffman)::RegHoffman
    r = RegHoffman()
    r.terms[1] = w
    return r
end
function RegHoffman(w::Word)::RegHoffman
    r = RegHoffman()
    r.terms[1] = HoffmanWordtoHoffman(w)
    return r
end

RegHoffman(r::RegHoffman)::RegHoffman = r
function RegHoffman(a::NN)::RegHoffman
    r = RegHoffman()
    r.terms[0] = Hoffman(a)
    return r
end
# TODO: RegHoffman for hrm shf mpl