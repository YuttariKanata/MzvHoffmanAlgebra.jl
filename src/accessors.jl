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

#function uppernumber: 数字を上付き文字に変換する
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

#show: Operatorの出力の見栄えをよくする
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
    if sym == :HoftoIdx
        return HoffmanWordtoIndex(w)
    elseif sym == :IdxtoHof
        return IndexWordtoHoffman(w)
    elseif sym == :HoftoHof
        return HoffmanWordtoHoffman(w)
    elseif sym == :IdxtoIdx
        return IndexWordtoIndex(w)
    elseif sym == :tovec
        return collect(Int64,w)
    else
        return getfield(w,sym)
    end
end

function getproperty(r::Poly, sym::Symbol)
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

# show (Word, Hoffman Index)

function show(io::IO, ::MIME"text/plain", w::Word)
    show(io, MIME("text/plain"), w.t)
end
function show(io::IO, ::MIME"text/plain", w::Hoffman)
    if iszero(w)
        print(io, "0")
        return
    end

    first = true
    lid = length(w.terms)
    for (k, (word, coeff)) in enumerate(w.terms)
        if k >= _OMIT_COUNTS && lid >= _OMIT_COUNTS + 100
            if k <= lid-30
                if k == _OMIT_COUNTS
                    printstyled("...[$(lid - _OMIT_COUNTS - 29) terms]...",color = 11)
                end
                continue
            end
        end
        # ±の表示
        if first
            first = false
            if coeff < Clong(0)
                print(io, "- ")
                coeff = -coeff
            end
            if coeff == Culong(1)
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
function show(io::IO, w::Hoffman)
    if iszero(w)
        print(io, "0")
        return
    end

    first = true
    for (word, coeff) in w.terms
        # ±の表示
        if first
            first = false
            if coeff < Clong(0)
                print(io, "- ")
                coeff = -coeff
            end
            if coeff == Culong(1)
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
function show(io::IO, ::MIME"text/plain", idx::Index)
    if iszero(idx)
        print(io, "0")
        return
    end

    first = true
    lid = length(idx.terms)
    for (k, (word, coeff)) in enumerate(idx.terms)
        if k >= _OMIT_COUNTS && lid >= _OMIT_COUNTS + 100
            if k <= lid-30
                if k == _OMIT_COUNTS
                    printstyled("...[$(lid- _OMIT_COUNTS - 29) terms]...",color = 11)
                end
                continue
            end
        end
        # ±の表示
        if first
            first = false
            if coeff < Culong(0)
                print(io, "- ")
                coeff = -coeff
            end
            if coeff == Culong(1)
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
function show(io::IO, idx::Index)
    if iszero(idx)
        print(io, "0")
        return
    end

    first = true
    for (word, coeff) in idx.terms
        # ±の表示
        if first
            first = false
            if coeff < Culong(0)
                print(io, "- ")
                coeff = -coeff
            end
            if coeff == Culong(1)
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

# sortedprint (Hoffman, Index)

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

# naturalshow (Hoffman, Index)

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
                if isempty(word)
                    print(io, coeff_to_str(-coeff), "")
                end
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
function naturalshow(io::IO, w::Index, f::Bool)
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
                if isempty(word)
                    print(io, coeff_to_str(-coeff), "")
                end
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

        if isempty(word)
            print(io, "")
        else
            print(io, "[", join(word,", "), "]" )
        end
    end
end
function naturalshow(io::IO, r::Poly{Rational{BigInt}})
    if iszero(r)
        print(io, "0")
        return
    end

    first = true
    degs = sort(collect(keys(r.terms)), rev=true)
    for d in degs
        coeff = r.terms[d]
        if first
            first = false
            if coeff == Clong(-1)
                print(io, "- ")
                if d == 0
                    print(io, coeff_to_str(-coeff), "")
                end
            elseif coeff == Clong(1)
                if d == 0
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
            if coeff != Culong(1) || d == 0
                print(io, coeff_to_str(coeff), "")
            end
        end

        if d > 1
            print(io, "T^", d)
        elseif d == 1
            print(io, "T")
        end
    end
end

# show, sortedprint (Poly)

function show(io::IO, ::MIME"text/plain", r::Poly{A}) where A
    if !(A === Hoffman || A === Index || A === Rational{BigInt})
        @warn "Poly{$A}に対するshowはまだ未実装です。"
        println(r)
        return
    end
    if iszero(r)
        print(io, "0")
        return
    end
    if A === Rational{BigInt}
        naturalshow(io, r)
        return
    end
    degs = sort(collect(keys(r.terms)), rev=true)

    first = true
    for d in degs
        ci = r.terms[d]
        if d == 0
            naturalshow(io, ci, first)
            continue
        else
            if is_monomial(ci)
                if !isone(ci)
                    naturalshow(io,ci,first)
                end
            else
                print(io," + (")
                naturalshow(io, ci, first)
                print(io,")")
            end
        end
        if d == 1
            print(io, "T")
        else
            print(io, "T^",d)
        end
        first = false
    end
end
show(io::IO, r::Poly) = show(io, MIME("text/plain"), r)
function sortedprint(r::Poly{A}) where A
    if !(A === Hoffman || A === Index)
        @warn "Poly{$A}に対するsortedshowはまだ未実装です。"
        println(r)
        return
    end
    if iszero(r)
        println("Poly{$A} with no entry:")
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
        println("Poly{$A} with 1 entry:")
    else
        println("Poly{$A} with $(length(r.terms)) entries:")
    end

    if A === Hoffman
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
    elseif A === Index
        for d in degs
            println("T^",d,":")
            ch = r.terms[d]
            # Hoffman の流用
            words = keys(ch.terms)
            
            sorted_word = sort(collect(words); by = t -> (length(t), t))
    
            for wo in sorted_word
                c = coeff_to_str(ch.terms[wo])
                println(lpad(c,number_space)," ","[", join(wo,", "), "]")
            end
        end
    end
end
