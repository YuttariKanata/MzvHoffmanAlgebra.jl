#[ operator.jl ]#

# Defines Operator-related functions and constants


#=
export left_act, right_act, ⬆️, ➡️, ⬇️, ⬅️, ➖, up, right, down, left, minus, τ, 𖼷, η, ⋁, φ, ∂
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
            key = word(op.cnt)
        else
            key = word(w[1]+op.cnt,w[2:end]...)
        end
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function left_act(op::OpDown, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = word(-op.cnt)
        else
            key = word(w[1]-op.cnt,w[2:end]...)
        end
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function left_act(op::OpLeft, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        key = word(1)^op.cnt * w
        r.terms[key] = get(r.terms, key, 0) + c
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
            key = word(op.cnt)
        else
            key = word(w[1:end-1]..., w[end]+op.cnt)
        end
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function right_act(op::OpDown, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        if isone(w)
            key = word(-op.cnt)
        else
            key = word(w[1:end-1]..., w[end]-op.cnt)
        end
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end
function right_act(op::OpRight, t::Index)::Index
    r = Index()
    for (w,c) in t.terms
        key = w * word(1)^op.cnt
        r.terms[key] = get(r.terms, key, 0) + c
    end
    return r
end

function *(a::Operator, b::Operator)::Operator
    r = copy(a)
    append_clean!(r,b)
    return r
end
function *(a::Operator, bn::AbstractOp)::Operator
    r = copy(a)
    rlast = (isempty(r.ops) ? nothing : r.ops[end])

    if bn isa OpDeriv
        b = OpDeriv(bn.n,bn.cnt)
    elseif bn isa OpDown
        b = OpUp(-bn.cnt)
    elseif bn isa OpTau
        b = OpTau(bn.cnt & 1)
    else
        b = typeof(bn)(bn.cnt)
    end

    tb = typeof(b)
    ta = typeof(rlast)
    if ta !== tb
        push!(r.ops,b)
    elseif ta === OpTau
        r.ops[end].cnt = xor(b.cnt, rlast.cnt) & 1
    elseif ta === OpDeriv
        if b.n == rlast.n
            r.ops[end].cnt += b.cnt
        end
    else
        r.ops[end] += b.cnt
    end

    if !isempty(r.ops) && r.ops[end].cnt == 0
        pop!(r.ops)
    end

    return r
end
function *(a::Operator,bn::Type{<:AbstractOp})::Operator
    r = copy(a)
    rlast = (isempty(r.ops) ? nothing : r.ops[end])

    if bn === OpDown
        b = OpUp
    end

    ta = typeof(rlast)
    if ta !== b
        push!(r.ops, one(b))
    end
end

function ^(op::Type{<:AbstractOp}, n::Integer)::Operator
    if n == 0
        return Operator()
    else
        return Operator((op)(n))
    end 
end
function ^(op::Operator, n::Integer)::Operator
    if n < 0
        return throw(DomainError(n, "negative power is not supported"))
    elseif n == 0
        return one(Operator)
    else
        return op * op^(n-1)
    end
end