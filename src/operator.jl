#[ operator.jl ]#

# Defines Operator-related functions and constants


#=
export left_act, right_act, ⬆️, ➡️, ⬇️, ⬅️, ➖, up, right, down, left, minus, τ, 𖼷, η, ⋁, φ, ∂,
       WordtoOperator
=#



"""
###################################################################################################
                                        Operators
###################################################################################################
"""


###################################################################################################
############## Operator Constants #################################################################


const ⬆️::DataType    = OpUp
const ➡️::DataType    = OpRight
const ⬇️::DataType    = OpDown
const ⬅️::DataType    = OpLeft
const ➖::DataType    = OpMinus
const up::DataType    = OpUp
const right::DataType = OpRight
const down::DataType  = OpDown
const left::DataType  = OpLeft
const minus::DataType = OpMinus
const τ::DataType     = OpTau
const 𖼷::DataType     = OpTau
const η::DataType     = OpEta
const ⋁::DataType     = OpEta
const φ::DataType     = OpPhi
const ∂::DataType     = OpDeriv


###################################################################################################
############## Operator function ##################################################################

############################## left_act ##############################
function left_act(op::OpMinus, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        key = w[(1+op.cnt):end]
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function left_act(op::OpUp, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = Word(op.cnt)
        else
            key = Word(w[1]+op.cnt,w[2:end]...)
        end
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function left_act(op::OpDown, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = Word(-op.cnt)
        else
            key = Word(w[1]-op.cnt,w[2:end]...)
        end
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function left_act(op::OpLeft, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        key = Word(1)^op.cnt * w
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function left_act(op::OpTau, t::Index)::Index
    r = copy(t)
    if op.cnt & 1 == 1
        return dual(r)
    else
        return r
    end
end
function left_act(op::OpEta, t::Index)::Index
    r = copy(t)
    if op.cnt == 0
        return r
    else
        for _ in 1:op.cnt
            r = Hoffman_dual(r)
        end
    end
    return r
end
function left_act(op::OpPhi, t::Index)::Index
    r = copy(t)
    if op.cnt == 0
        return r
    else
        for _ in 1:op.cnt
            r = Landen_dual(r)
        end
    end
    return r
end

left_act(op::OpMinus, w::Word)::Word = w[1+op.cnt:end]
left_act(op::OpUp, w::Word)::Word    = isone(w) ? Word(op.cnt)  : Word(w[1]+op.cnt,w[2:end]...)
left_act(op::OpDown, w::Word)::Word  = isone(w) ? Word(-op.cnt) : Word(w[1]-op.cnt,w[2:end]...)
left_act(op::OpLeft, w::Word)::Word  = Word(1)^op.cnt * w
left_act(op::OpTau, w::Word)::Word   = op.cnt & 1 == 1 ? monomial_dual_i(w) : w
function left_act(op::OpEta, w::Word)::Word
    if op.cnt == 0
        return w
    end
    r = w
    for _ in 1:op.cnt
        r = monomial_hof_dual_i(r)
    end
    return r
end
function left_act(op::OpPhi, w::Word)::Index
    if op.cnt == 0
        return w
    end
    r = w
    for _ in 1:op.cnt
        r = Landen_dual(r)
    end
    return r
end

############################## right_act ##############################
function right_act(op::OpMinus, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        key = w[1:end-op.cnt]
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function right_act(op::OpUp, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = Word(op.cnt)
        else
            key = Word(w[1:end-1]..., w[end]+op.cnt)
        end
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function right_act(op::OpDown, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = Word(-op.cnt)
        else
            key = Word(w[1:end-1]..., w[end]-op.cnt)
        end
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function right_act(op::OpRight, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        key = w * Word(1)^op.cnt
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end

right_act(op::OpMinus, w::Word)::Word = w[1:end-op.cnt]
right_act(op::OpUp, w::Word)::Word    = isone(w) ? Word(op.cnt)  : Word(w[1:end-1]..., w[end]+op.cnt)
right_act(op::OpDown, w::Word)::Word  = isone(w) ? Word(-op.cnt) : Word(w[1:end-1]..., w[end]-op.cnt)
right_act(op::OpRight, w::Word)::Word = w * Word(1)^op.cnt

############################## multiplication ##############################
function *(a::Operator, b::Operator)::Operator
    r = copy(a)
    append_clean!(r,b)
    return r
end
function *(a::Operator, bn::AbstractOp)::Operator
    r = copy(a)
    rlast = (isempty(r.ops) ? nothing : r.ops[end])
    b = copy(bn)

    tb = typeof(b)
    ta = typeof(rlast)
    if ta !== tb
        push!(r.ops,b)
    elseif ta === OpTau
        r.ops[end].cnt = xor(b.cnt, rlast.cnt) & 1
    elseif ta === OpDeriv
        if b.n == rlast.n
            r.ops[end].cnt += b.cnt
        else
            push!(r.ops,b)
        end
    else
        r.ops[end].cnt += b.cnt
    end

    if !isempty(r.ops) && r.ops[end].cnt == 0
        pop!(r.ops)
    end

    return r
end
function *(an::AbstractOp, b::Operator)::Operator
    r = copy(b)
    rfirst = (isempty(r.ops) ? nothing : r.ops[1])
    a = copy(an)

    ta = typeof(a)
    tb = typeof(rfirst)
    if ta !== tb
        pushfirst!(r.ops,a)
    elseif ta === OpTau
        r.ops[1].cnt = xor(a.cnt, rfirst.cnt) & 1
    elseif ta === OpDeriv
        if a.n == rfirst.n
            r.ops[1].cnt += a.cnt
        else
            pushfirst!(r.ops,a)
        end
    else
        r.ops[1].cnt += a.cnt
    end

    if !isempty(r.ops) && r.ops[1].cnt == 0
        popfirst!(r.ops)
    end

    return r
end
# 一個一個いちいちめんどくさいのでやめる
*(a::Operator,bn::Type{<:AbstractOp})::Operator            = *(a,(bn)(1))
*(an::Type{<:AbstractOp},b::Operator)::Operator            = *((an)(1),b)
*(an::AbstractOp,bn::AbstractOp)::Operator                 = *(Operator(an),bn)
*(an::AbstractOp,bn::Type{<:AbstractOp})::Operator         = *(Operator(an),bn)
*(an::Type{<:AbstractOp},bn::AbstractOp)::Operator         = *(an,Operator(bn))
*(an::Type{<:AbstractOp},bn::Type{<:AbstractOp})::Operator = *(Operator(an),(bn)(1))

############################## powers ##############################
function ^(op::Type{<:AbstractOp}, n::Integer)::Operator
    if n == 0
        return Operator()
    else
        return Operator((op)(n))
    end 
end
function ^(op::AbstractOp, m::Integer)::Operator
    r = Operator()
    t = typeof(op)
    if t === OpDeriv
        push!(r.ops,OpDeriv(op.n,op.cnt*m))
    elseif t === OpDown
        push!(r.ops,OpUp(-op.cnt*m))
    else
        push!(r.ops,t(op.cnt*m))
    end
    return r
end
function ^(op::Operator, n::Integer)::Operator
    if n < 0
        return throw(DomainError(n, "negative power is not supported"))
    elseif n == 0
        return one(Operator)
    elseif n == 1
        return op
    else
        return op * op^(n-1)
    end
end

###################################################################################################
############## Operator functions related to word #################################################
function WordtoOperator(w::Word)::Operator
    r = Operator()
    lw = lastindex(w)
    if lw == 0
        return r
    end
    v = Vector{AbstractOp}(undef,lw*2)
    for i in 1:lw
        v[2*i-1] = OpRight(1)
        v[2*i] = OpUp(w[i]-1)
    end
    append_clean!(r,v)
    return r
end

###################################################################################################
############## Action function with Index #########################################################
function *(op::Operator, idx::Index)::Index
    ridx = copy(idx)
    if isone(op)
        return idx
    end
    for i in lastindex(op):-1:1
        ridx = left_act(op.ops[i],ridx)
    end
    return ridx
end
*(op::AbstractOp, idx::Index)::Index = left_act(op,copy(idx))
*(op::Type{<:AbstractOp}, idx::Index)::Index = left_act((op)(1),copy(idx))
function *(idx::Index, op::Operator)::Index
    ridx = copy(idx)
    if isone(op)
        return idx
    end
    for i in 1:lastindex(op)
        ridx = right_act(op.ops[i],ridx)
    end
    return ridx
end
*(idx::Index, op::AbstractOp)::Index = right_act(op,copy(idx))
*(idx::Index, op::Type{<:AbstractOp})::Index = right_act((op)(1),copy(idx))

function *(op::Operator, w::Word)::Word
    rw = w
    if isone(op)
        return rw
    end
    for i in lastindex(op):-1:1
        rw = left_act(op.ops[i],rw)
    end
    return rw
end
*(op::AbstractOp, w::Word)::Word = left_act(op,w)
*(op::Type{<:AbstractOp}, w::Word)::Word = left_act((op)(1),w)
function *(w::Word, op::Operator)::Word
    rw = w
    if isone(op)
        return w
    end
    for i in 1:lastindex(op)
        rw = right_act(op.ops[i],rw)
    end
    return rw
end
*(w::Word, op::AbstractOp)::Word = right_act(op,w)
*(w::Word, op::Type{<:AbstractOp})::Word = right_act((op)(1),w)