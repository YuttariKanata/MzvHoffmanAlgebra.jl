#[ operators.jl ]#

# Defines Operator-related functions and constants

import Base: *, ^

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


const ⬆️    = OpUp
const ➡️    = OpRight
const ⬇️    = OpDown
const ⬅️    = OpLeft
const ➖    = OpMinus
const up    = OpUp
const right = OpRight
const down  = OpDown
const left  = OpLeft
const minus = OpMinus
const τ     = OpTau
const 𖼷     = OpTau
const η     = OpEta
const ⋁     = OpEta
const φ     = OpPhi
const ∂     = OpDeriv


###################################################################################################
############## Operator function ##################################################################

############################## left_act ##############################

function left_act(op::OpMinus, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        # w is IndexWord, w.t is Tuple
        if length(w) > op.cnt
            key = IndexWord(w.t[(1+op.cnt):end])
            _add!(r.terms, key, c)
        end
    end
    return r
end

function left_act(op::OpUp, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = IndexWord((op.cnt,))
        else
            # Add op.cnt to the first index
            key = IndexWord((w.t[1]+op.cnt, w.t[2:end]...))
        end
        _add!(r.terms, key, c)
    end
    return r
end

function left_act(op::OpDown, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = IndexWord((-op.cnt,))
        else
            key = IndexWord((w.t[1]-op.cnt, w.t[2:end]...))
        end
        _add!(r.terms, key, c)
    end
    return r
end

function left_act(op::OpLeft, t::Index)::Index
    r = Index()
    prefix = ntuple(i -> 1, op.cnt)
    for (w,c) in t.terms
        key = IndexWord((prefix..., w.t...))
        _add!(r.terms, key, c)
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
        r = Hoffman(t)
        for _ in 1:op.cnt
            r = Hoffman_dual(r)
        end
    end
    return Index(r)
end

function left_act(op::OpPhi, t::Index)::Index
    r = copy(t)
    if op.cnt == 0
        return r
    else
        r = Hoffman(t)
        for _ in 1:op.cnt
            r = Landen_dual(r)
        end
    end
    return Index(r)
end

function left_act(op::OpDeriv, t::Index)::Index
    r = copy(t)
    if op.cnt == 0
        return r
    else
        r = Hoffman(t)
        for _ in 1:op.cnt
            r = dell(r,op.n)
        end
    end
    return Index(r)
end

# --- Action on IndexWord ---

left_act(op::OpMinus, w::IndexWord)::IndexWord = length(w) > op.cnt ? IndexWord(w.t[1+op.cnt:end]) : IndexWord()
left_act(op::OpUp, w::IndexWord)::IndexWord    = isone(w) ? IndexWord((op.cnt,)) : IndexWord((w.t[1]+op.cnt, w.t[2:end]...))
left_act(op::OpDown, w::IndexWord)::IndexWord  = isone(w) ? IndexWord((-op.cnt,)) : IndexWord((w.t[1]-op.cnt, w.t[2:end]...))
left_act(op::OpLeft, w::IndexWord)::IndexWord  = IndexWord((ntuple(i->1, op.cnt)..., w.t...))

function left_act(op::OpTau, w::IndexWord)::IndexWord
    if op.cnt & 1 == 1
        monomial_dual(w)
    else
        return w
    end
end

function left_act(op::OpEta, w::IndexWord)::IndexWord
    if op.cnt == 0
        return w
    end
    r = HoffmanWord(w)
    for _ in 1:op.cnt
        r = monomial_hof_dual(r)
    end
    return IndexWord(r)
end

function left_act(op::OpPhi, w::IndexWord)::Index
    if op.cnt == 0
        return Index(w)
    end
    r = Hoffman(w)
    for _ in 1:op.cnt
        r = Landen_dual(r)
    end
    return Index(r)
end

function left_act(op::OpDeriv, w::IndexWord)::Index
    if op.cnt == 0
        return Index(w)
    end
    r = Hoffman(w)
    for _ in 1:op.cnt
        r = dell(r,op.n)
    end
    return Index(r)
end

############################## right_act ##############################

function right_act(op::OpMinus, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if length(w) > op.cnt
            key = IndexWord(w.t[1:end-op.cnt])
            _add!(r.terms, key, c)
        end
    end
    return r
end

function right_act(op::OpUp, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = IndexWord((op.cnt,))
        else
            key = IndexWord((w.t[1:end-1]..., w.t[end]+op.cnt))
        end
        _add!(r.terms, key, c)
    end
    return r
end

function right_act(op::OpDown, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = IndexWord((-op.cnt,))
        else
            key = IndexWord((w.t[1:end-1]..., w.t[end]-op.cnt))
        end
        _add!(r.terms, key, c)
    end
    return r
end

function right_act(op::OpRight, t::Index)::Index
    r = Index()
    suffix = ntuple(i -> 1, op.cnt)
    for (w,c) in t.terms
        key = IndexWord((w.t..., suffix...))
        _add!(r.terms, key, c)
    end
    return r
end

# --- Action on IndexWord ---

right_act(op::OpMinus, w::IndexWord)::IndexWord = length(w) > op.cnt ? IndexWord(w.t[1:end-op.cnt]) : IndexWord()
right_act(op::OpUp, w::IndexWord)::IndexWord    = isone(w) ? IndexWord((op.cnt,)) : IndexWord((w.t[1:end-1]..., w.t[end]+op.cnt))
right_act(op::OpDown, w::IndexWord)::IndexWord  = isone(w) ? IndexWord((-op.cnt,)) : IndexWord((w.t[1:end-1]..., w.t[end]-op.cnt))
right_act(op::OpRight, w::IndexWord)::IndexWord = IndexWord((w.t..., ntuple(i->1, op.cnt)...))

############################## multiplication ##############################

function *(a::Operator, b::Operator)::Operator
    # Operator constructor calls clean, so we just concatenate
    return Operator([a.ops; b.ops])
end

*(a::Operator, bn::AbstractOp)::Operator = Operator([a.ops; bn])
*(an::AbstractOp, b::Operator)::Operator = Operator([an; b.ops])

*(a::Operator, bn::Type{<:AbstractOp})::Operator            = a * (bn)(1)
*(an::Type{<:AbstractOp}, b::Operator)::Operator            = (an)(1) * b
*(an::AbstractOp, bn::AbstractOp)::Operator                 = Operator([an, bn])
*(an::AbstractOp, bn::Type{<:AbstractOp})::Operator         = Operator([an, (bn)(1)])
*(an::Type{<:AbstractOp}, bn::AbstractOp)::Operator         = Operator([(an)(1), bn])
*(an::Type{<:AbstractOp}, bn::Type{<:AbstractOp})::Operator = Operator([(an)(1), (bn)(1)])

############################## powers ##############################

function ^(op::Type{<:AbstractOp}, n::Integer)::Operator
    if n == 0
        return Operator()
    else
        return Operator([(op)(n)])
    end 
end

function ^(op::AbstractOp, m::Integer)::Operator
    # For single abstract op, we can often simplify directly
    t = typeof(op)
    if t === OpDeriv
        return Operator([OpDeriv(op.n, op.cnt * m)])
    elseif t === OpDown
        return Operator([OpUp(-op.cnt * m)])
    else
        return Operator([(t)(op.cnt * m)])
    end
end

function ^(op::Operator, n::Integer)::Operator
    if n < 0
        throw(DomainError(n, "negative power is not supported"))
    elseif n == 0
        return one(Operator)
    elseif n == 1
        return op
    else
        # Binary exponentiation could be used, but simple recursion is fine for small n
        return op * op^(n-1)
    end
end

###################################################################################################
############## Operator functions related to word #################################################

function WordtoOperator(w::IndexWord)::Operator
    r = Operator()
    lw = length(w)
    if lw == 0
        return r
    end
    
    # Construct the operator sequence
    # w = (k1, k2, ...) -> OpRight(1) * OpUp(k1-1) * OpRight(1) * OpUp(k2-1) ...
    ops = Vector{AbstractOp}(undef, lw * 2)
    for i in 1:lw
        ops[2*i-1] = OpRight(1)
        ops[2*i]   = OpUp(w[i]-1)
    end
    
    # Use clean to normalize
    return Operator(ops)
end

###################################################################################################
############## Action function with Index #########################################################

function *(op::Operator, idx::Index)::Index
    ridx = copy(idx)
    if isone(op)
        return ridx
    end
    # Apply operators from right to left (last to first in list)
    for i in lastindex(op.ops):-1:1
        ridx = left_act(op.ops[i], ridx)
    end
    return ridx
end

*(op::AbstractOp, idx::Index)::Index = left_act(op, copy(idx))
*(op::Type{<:AbstractOp}, idx::Index)::Index = left_act((op)(1), copy(idx))

function *(idx::Index, op::Operator)::Index
    ridx = copy(idx)
    if isone(op)
        return ridx
    end
    # Apply operators from left to right
    for i in 1:lastindex(op.ops)
        ridx = right_act(op.ops[i], ridx)
    end
    return ridx
end

*(idx::Index, op::AbstractOp)::Index = right_act(op, copy(idx))
*(idx::Index, op::Type{<:AbstractOp})::Index = right_act((op)(1), copy(idx))

# --- Action on Word (returning Word or Index) ---

function *(op::Operator, w::IndexWord)::Union{IndexWord, Index}
    res = w
    if isone(op)
        return res
    end
    
    for i in lastindex(op.ops):-1:1
        res = left_act(op.ops[i], res)
    end
    return res
end

*(op::AbstractOp, w::IndexWord)::Union{IndexWord, Index} = left_act(op, w)
*(op::Type{<:AbstractOp}, w::IndexWord)::Union{IndexWord, Index} = left_act((op)(1), w)

function *(w::IndexWord, op::Operator)::IndexWord
    # Right actions usually preserve the single word structure for standard operators
    rw = w
    if isone(op)
        return w
    end
    for i in 1:lastindex(op.ops)
        rw = right_act(op.ops[i], rw)
    end
    return rw
end

*(w::IndexWord, op::AbstractOp)::IndexWord = right_act(op, w)
*(w::IndexWord, op::Type{<:AbstractOp})::IndexWord = right_act((op)(1), w)