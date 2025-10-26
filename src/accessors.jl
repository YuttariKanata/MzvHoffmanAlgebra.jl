#[ accessors.jl ]#

# This file defines functions that extend Base-derived functions


import Base: show, getproperty

#=
export upper_represent, sortedprint
=#

"""
###################################################################################################
                                        Operators
###################################################################################################
"""


###################################################################################################
############## about representation ###############################################################
const _UpperNumber_Table = Dict{Int,Char}(
    0 => '⁰',
    1 => '¹',
    2 => '²',
    3 => '³',
    4 => '⁴',
    5 => '⁵',
    6 => '⁶',
    7 => '⁷',
    8 => '⁸',
    9 => '⁹',
    -1 => '⁻',
)

@inline function uppernumber(n::Integer)::String
    ds = digits(abs(n))
    l = lastindex(ds)
    v = Vector{Char}(undef,l)
    for i in 1:l
        v[i] = _UpperNumber_Table[ds[l+1-i]]
    end
    if signbit(n)
        pushfirst!(v,_UpperNumber_Table[-1])
    end
    return join(v)
end

const _upper_represent = Base.RefValue{Bool}(false)
upper_represent(val::Bool) = (_upper_represent[] = val)

function show(io::IO, ::MIME"text/plain", op::Operator)
    if isempty(op.ops)
        print(io, "∅")
        return
    end
    opl = lastindex(op.ops)
    parts = Vector{String}(undef,opl)
    if _upper_represent[]
        for i in 1:opl
            ot = typeof(op.ops[i])
            k = op.ops[i].cnt
            if ot === OpDeriv
                if k == 1
                    parts[i] = "∂" * string(op.ops[i].n)
                else
                    parts[i] = "∂" * string(op.ops[i].n) * uppernumber(k)
                end
            else
                if k == 1
                    parts[i] =  _Operator_to_String_Table[ot]
                else
                    parts[i] =  _Operator_to_String_Table[ot] * uppernumber(op.ops[i].cnt)
                end
            end
        end
    else
        for i in 1:opl
            ot = typeof(op.ops[i])
            k = op.ops[i].cnt
            if ot === OpDeriv
                if k == 1
                    parts[i] = "∂" * string(op.ops[i].n)
                elseif k > 0
                    parts[i] = "∂" * string(op.ops[i].n) * "^" * string(k)
                else
                    parts[i] = "∂" * string(op.ops[i].n) * "^(" * string(k) * ")"
                end
            else
                if k == 1
                    parts[i] =  _Operator_to_String_Table[ot]
                elseif k > 0
                    parts[i] =  _Operator_to_String_Table[ot] * "^" * string(op.ops[i].cnt)
                else
                    parts[i] =  _Operator_to_String_Table[ot] * "^(" * string(op.ops[i].cnt) * ")"
                end
            end
        end
    end

    print(io, join(parts, " * "))
end
show(io::IO, op::Operator) = show(io, MIME("text/plain"), op)
function show(io::IO, ::MIME"text/plain", op::AbstractOp)
    if op isa OpDeriv
        if _upper_represent[]
            print(io,"∂",op.n,uppernumber(op.cnt))
        else
            print(io,"∂",op.n,"^",op.cnt)
        end
    else
        if _upper_represent[]
            print(io,_Operator_to_String_Table[typeof(op)],uppernumber(op.cnt))
        else
            print(io,_Operator_to_String_Table[typeof(op)],"^",op.cnt)
        end
    end
end
show(io::IO, op::AbstractOp) = show(io, MIME("text/plain"), op)









"""
###################################################################################################
                                        Hoffman MZV
###################################################################################################
"""


###################################################################################################
############## get property #######################################################################


function getproperty(w::Union{MonoIndex,Hoffman,Index}, sym::Symbol)
    if sym == :toIndex
        return Index(w)
    elseif sym == :toHoffman
        return Hoffman(w)
    elseif sym == :toMonoIndex
        return MonoIndex(w)
    elseif sym == :toWord
        return Word(w)
    elseif sym == :sortshow
        sortedprint(w)
    else
        return getfield(w,sym)
    end
end

function getproperty(w::Word, sym::Symbol)
    if sym == :toindex
        return idxdprs(collect(w))
    elseif sym == :toword
        if get_index_orientation()
            return Word(idxprs(w))
        else
            return Word(idxprs_r(w))
        end
    elseif sym == :tovec
        return collect(Int64,w)
    elseif sym == :last
        return last_letter(w)
    elseif sym == :rmlast
        return remove_lastletter(w)
    else
        return getfield(w,sym)
    end
end

function getproperty(r::RegHoffman, sym::Symbol)
    if sym == :sortshow
        sortedprint(r)
    else
        return getfield(r,sym)
    end
end

###################################################################################################
############## about representation ###############################################################

@inline function coeff_to_str(r::Rational{BigInt})::String
    if r.den == 1
        return string(r.num)
    end
    return string(r.num) * "/" * string(r.den)
end
@inline function word_to_str(w::Word)::String
    s = ""
    for i in w
        if i == 1
            s *= "x"
        elseif i == 2
            s *= "y"
        end
    end
    return s
end
function show(io::IO, ::MIME"text/plain", w::Word)
    show(io, MIME("text/plain"), w.t)
end
function show(io::IO, ::MIME"text/plain", w::Hoffman)
    if iszero(w)
        print(io, "0")
        return
    end

    first = true
    for (word, coeff) in w.terms
        # ±の表示
        if first
            first = false
            if coeff == Clong(-1)
                print(io, "- ")
            elseif coeff == Culong(1)
                if isempty(word)
                    print(io, coeff_to_str(coeff), "")
                end
            else
                print(io, coeff_to_str(coeff), "")
            end
        else
            if coeff < 0
                print(io, " - ")
                coeff = -coeff
            else
                print(io, " + ")
            end
            if coeff != Culong(1) || isempty(word)
                print(io, coeff_to_str(coeff), "")
            end
        end

        print(io, word_to_str(word) )
    end
end
show(io::IO, w::Hoffman) = show(io, MIME("text/plain"), w)
function show(io::IO, ::MIME"text/plain", idx::Index)
    if iszero(idx)
        print(io, "0")
        return
    end

    first = true
    for (word, coeff) in idx.terms
        # ±の表示
        if first
            first = false
            if coeff == Clong(-1)
                print(io, "- ")
            elseif coeff == Culong(1)
                if isempty(word)
                    print(io, coeff_to_str(coeff), "")
                end
            else
                print(io, coeff_to_str(coeff), "")
            end
        else
            if coeff < 0
                print(io, " - ")
                coeff = -coeff
            else
                print(io, " + ")
            end
            if coeff != 1//1 || isempty(word)
                print(io, coeff_to_str(coeff), "")
            end
        end

        if !isempty(word)
            print(io, "[", join(word,", "), "]" )
        end
    end
end
show(io::IO, idx::Index) = show(io, MIME("text/plain"), idx)

function sortedprint(w::Hoffman)
    if iszero(w)
        println("Hoffman with no entry:")
        println("    0")
        return
    end

    words = keys(w.terms)
    coeffs = values(w.terms)
    
    max_coeff = maximum(t-> ndigits(t.den)+ndigits(t.num),coeffs)
    number_space = max_coeff+3
    
    sorted_word = sort(collect(words); by = t -> (length(t), t))

    if length(w.terms) == 1
        println("Hoffman with $(length(w.terms)) entry:")
    else
        println("Hoffman with $(length(w.terms)) entries:")
    end
    for wo in sorted_word
        c = coeff_to_str(w.terms[wo])
        println(lpad(c,number_space)," ",word_to_str(wo))
    end
end
function sortedprint(idx::Index)
    if iszero(idx)
        println("Index with no entry:")
        println("    0 []")
        return
    end

    words = keys(idx.terms)
    coeffs = values(idx.terms)
    
    max_coeff = maximum(t-> ndigits(t.den)+ndigits(t.num),coeffs)
    number_space = max_coeff+3
    
    sorted_word = sort(collect(words); by = t -> (sum(t), t))

    if length(idx.terms) == 1
        println("Index with $(length(idx.terms)) entry:")
    else
        println("Index with $(length(idx.terms)) entries:")
    end
    for wo in sorted_word
        c = coeff_to_str(idx.terms[wo])
        println(lpad(c,number_space)," [",join(wo,", "),"]")
    end
end
function naturalshow(io::IO, w::Hoffman, f::Bool)
    if iszero(w)
        print(io, "0")
        return
    end

    first = f
    for (word, coeff) in w.terms
        # ±の表示
        if first
            first = false
            if coeff == Clong(-1)
                print(io, "- ")
            elseif coeff == Culong(1)
                if isempty(word)
                    print(io, coeff_to_str(coeff), "")
                end
            else
                print(io, coeff_to_str(coeff), "")
            end
        else
            if coeff < 0
                print(io, " - ")
                coeff = -coeff
            else
                print(io, " + ")
            end
            if coeff != Culong(1) || isempty(word)
                print(io, coeff_to_str(coeff), "")
            end
        end

        print(io, word_to_str(word) )
    end
end
function show(io::IO, ::MIME"text/plain", r::RegHoffman)
    if iszero(r)
        print(io, "0")
        return
    end
    degs = sort(collect(keys(r.terms)), rev=true)

    Sg = ""
    first = true
    for d in degs
        ch = r.terms[d]
        if d == 0
            print(io, Sg)
            show(io, MIME("text/plain"), ch)
            Sg = " + "
            continue
        else
            if is_monomial(ch)
                if !isone(ch)
                    naturalshow(io,ch,first)
                end
            else
                print(io, Sg, "(")
                show(io, MIME("text/plain"), ch)
                print(io, ")")
            end
        end
        Sg = " + "
        if d == 1
            print(io, "T")
        else
            print(io, "T^",d)
        end
        first = false
    end
end
show(io::IO, r::RegHoffman) = show(io, MIME("text/plain"), r)
function sortedprint(r::RegHoffman)
    if iszero(r)
        println("RegHoffman with no entry:")
        println("    0")
        return
    end

    number_space = 0
    degs = sort(collect(keys(r.terms)),rev=true)
    for (d,ch) in r.terms
        coeffs = values(ch.terms)
        max_coeff = maximum(t-> ndigits(t.den)+ndigits(t.num),coeffs)
        number_space = max(max_coeff+3,number_space)
    end

    if length(r.terms) == 1
        println("RegHoffman with 1 entry:")
    else
        println("RegHoffman with $(length(r.terms)) entries:")
    end

    for d in degs
        println("T^",d,":")
        ch = r.terms[d]
        # Hoffman の流用
        words = keys(ch.terms)
        
        sorted_word = sort(collect(words); by = t -> (length(t), t))

        for wo in sorted_word
            c = coeff_to_str(ch.terms[wo])
            println(lpad(c,number_space)," ",word_to_str(wo))
        end
    end
end