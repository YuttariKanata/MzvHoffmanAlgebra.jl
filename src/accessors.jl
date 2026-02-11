#[ accessors.jl ]#

# This file defines functions for string representation and visualization.

import Base: show, getproperty

"""
###################################################################################################
                                        Representation Helpers
###################################################################################################
"""

# Superscript mapping
const _UpperNumber_Table = Dict{Int,Char}(
    0 => '⁰', 1 => '¹', 2 => '²', 3 => '³', 4 => '⁴',
    5 => '⁵', 6 => '⁶', 7 => '⁷', 8 => '⁸', 9 => '⁹',
    -1 => '⁻'
)

@inline function uppernumber(n::Integer)::String
    ds = digits(abs(n))
    l = lastindex(ds)
    v = Vector{Char}(undef, l)
    for i in 1:l
        v[i] = _UpperNumber_Table[ds[l+1-i]]
    end
    if signbit(n)
        pushfirst!(v, _UpperNumber_Table[-1])
    end
    return join(v)
end

# Global toggle for upper case representation in Operators
const _upper_represent = Base.RefValue{Bool}(false)
"""
    upper_represent(val::Bool)

Sets whether to use superscript notation (e.g., ∂₁²) for operators in output.
"""
upper_represent(val::Bool) = (_upper_represent[] = val)

@inline function coeff_to_str(r::Rational{BigInt})::String
    if r.den == 1
        return string(r.num)
    end
    return string(r.num) * "/" * string(r.den)
end

@inline function word_to_str(w::HoffmanWord)::String
    s = ""
    for i in w.t
        if i == 0
            s *= "x"
        elseif i == 1
            s *= "y"
        else
            s *= "?" # Should not happen for valid HoffmanWord
        end
    end
    return s
end

@inline function word_to_str(w::IndexWord)::String
    isempty(w) && return ""
    return "Index(" * join(w.t, ",") * ")"
end

"""
###################################################################################################
                                        Show Methods
###################################################################################################
"""

# --- Operator ---

function show(io::IO, ::MIME"text/plain", op::Operator)
    if isempty(op.ops)
        print(io, "id") # Identity operator
        return
    end
    opl = lastindex(op.ops)
    parts = Vector{String}(undef, opl)
    
    use_upper = _upper_represent[]
    
    for i in 1:opl
        o = op.ops[i]
        ot = typeof(o)
        k = o.cnt
        
        base_str = ""
        if ot === OpDeriv
            base_str = "∂" * string(o.n)
        else
            base_str = get(_Operator_to_String_Table, ot, "?")
        end
        
        if use_upper
            if k == 1
                parts[i] = base_str
            else
                parts[i] = base_str * uppernumber(k)
            end
        else
            if k == 1
                parts[i] = base_str
            elseif k > 0
                parts[i] = base_str * "^" * string(k)
            else
                parts[i] = base_str * "^(" * string(k) * ")"
            end
        end
    end

    print(io, join(parts, " "))
end
show(io::IO, op::Operator) = show(io, MIME("text/plain"), op)

function show(io::IO, ::MIME"text/plain", op::AbstractOp)
    show(io, MIME("text/plain"), Operator([op]))
end
show(io::IO, op::AbstractOp) = show(io, MIME("text/plain"), op)

# --- Words ---

show(io::IO, w::HoffmanWord) = print(io, isempty(w) ? "1" : word_to_str(w))
show(io::IO, w::IndexWord) = print(io, isempty(w) ? "1" : word_to_str(w))

# --- Hoffman / Index (Linear Combinations) ---

function _show_linear_combination(io::IO, obj::Union{Hoffman, Index})
    if iszero(obj)
        print(io, "0")
        return
    end

    first = true
    lid = length(obj.terms)
    
    # Note: Iteration order of Dict is undefined. 
    # For deterministic output in tests/docs, sortedprint is preferred.
    # Here we just iterate for performance.
    
    for (k, (word, coeff)) in enumerate(obj.terms)
        # Omit if too long
        if k >= _OMIT_COUNTS && lid >= _OMIT_COUNTS + 100
            if k <= lid - 30
                if k == _OMIT_COUNTS
                    printstyled(io, " ...[$(lid - _OMIT_COUNTS - 29) terms]... ", color=:light_black)
                end
                continue
            end
        end

        # Sign handling
        if first
            first = false
            if coeff < 0
                print(io, "-")
                coeff = -coeff
            end
        else
            if coeff < 0
                print(io, " - ")
                coeff = -coeff
            else
                print(io, " + ")
            end
        end

        # Coefficient and Word
        is_unity = (coeff == 1)
        word_s = word_to_str(word)
        is_empty_word = isempty(word)

        if is_unity && !is_empty_word
            print(io, word_s)
        elseif is_empty_word
            print(io, coeff_to_str(coeff))
        else
            print(io, coeff_to_str(coeff), is_unity ? "" : " ", word_s)
        end
    end
end


show(io::IO, ::MIME"text/plain", h::Hoffman) = _show_linear_combination(io, h)
show(io::IO, h::Hoffman) = _show_linear_combination(io, h)

show(io::IO, ::MIME"text/plain", idx::Index) = _show_linear_combination(io, idx)
show(io::IO, idx::Index) = _show_linear_combination(io, idx)

# --- Poly ---

function naturalshow(io::IO, w::Union{Hoffman, Index}, first_term::Bool)
    if iszero(w)
        print(io, "0")
        return
    end

    first = first_term
    for (word, coeff) in w.terms
        if first
            first = false
            if coeff == -1
                print(io, "-")
                if isempty(word) print(io, "1") end
            elseif coeff == 1
                if isempty(word) print(io, "1") end
            else
                print(io, coeff_to_str(coeff))
            end
        else
            if coeff < 0
                print(io, " - ")
                coeff = -coeff
            else
                print(io, " + ")
            end
            if coeff != 1 || isempty(word)
                print(io, coeff_to_str(coeff))
            end
        end
        
        if !isempty(word)
            if coeff != 1 && coeff != -1
                print(io, " ")
            end
            print(io, word_to_str(word))
        end
    end
end

function show_poly_rat(io::IO, r::Poly{Rational{BigInt}})
    if iszero(r)
        print(io, "0")
        return
    end

    degs = sort(collect(keys(r.terms)), rev=true)
    degs = filter(d -> !iszero(r.terms[d]), degs)
    lid = length(degs)

    first = true
    for (k, d) in enumerate(degs)
        if k >= _OMIT_COUNTS && lid >= _OMIT_COUNTS + 100
            if k <= lid - 30
                if k == _OMIT_COUNTS
                    printstyled(io, " ...[$(lid - _OMIT_COUNTS - 29) terms]... ", color=:light_black)
                end
                continue
            end
        end

        coeff = r.terms[d]

        # Sign and separator
        if first
            first = false
            if coeff < 0
                print(io, "-")
                coeff = -coeff
            end
        else
            if coeff < 0
                print(io, " - ")
                coeff = -coeff
            else
                print(io, " + ")
            end
        end

        # Coefficient
        if coeff != 1 || d == 0
            print(io, coeff_to_str(coeff))
            if d > 0 print(io, " ") end
        end

        # Variable T
        if d == 1; print(io, "T"); elseif d > 1; print(io, "T^", d); end
    end
end

function show(io::IO, ::MIME"text/plain", r::Poly{A}) where A
    if iszero(r)
        print(io, "0")
        return
    end

    # Sort by degree descending
    if A === Rational{BigInt}
        show_poly_rat(io, r)
        return
    end

    degs = sort(collect(keys(r.terms)), rev=true)
    degs = filter(d -> !iszero(r.terms[d]), degs)
    lid = length(degs)

    first = true
    for (k, d) in enumerate(degs)
        if k >= _OMIT_COUNTS && lid >= _OMIT_COUNTS + 100
            if k <= lid - 30
                if k == _OMIT_COUNTS
                    printstyled(io, " ...[$(lid - _OMIT_COUNTS - 29) terms]... ", color=:light_black)
                end
                continue
            end
        end

        coeff = r.terms[d]
        
        if d == 0
            naturalshow(io, coeff, first)
            first = false
            continue
        end

        # For T^d terms
        is_mono = is_monomial(coeff)
        
        if isone(coeff)
            if !first print(io, " + ") end
        elseif is_mono && isone(-coeff)
            if first print(io, "-") else print(io, " - ") end
        elseif is_mono
            naturalshow(io, coeff, first)
            print(io, " ")
        else # not monomial
            if !first print(io, " + ") end
            print(io, "(")
            naturalshow(io, coeff, true)
            print(io, ")")
            print(io, " ")
        end

        if d == 1
            print(io, "T")
        elseif d > 1
            print(io, "T^", d)
        else
            # constant term was handled before
        end
        first = false
    end
end
show(io::IO, r::Poly) = show(io, MIME("text/plain"), r)

"""
    sortedprint(obj)

Prints the element with terms sorted by word length/weight for easier reading.
"""
function sortedprint(w::Hoffman)
    if iszero(w)
        println("Hoffman (0)")
        return
    end

    coeffs = values(w.terms)
    max_len = isempty(coeffs) ? 0 : maximum(c -> length(coeff_to_str(c)), coeffs)
    pad = max_len + 2
    
    # Sort by (length, word)
    sorted_keys = sort(collect(keys(w.terms)); by = t -> (length(t), t.t))
    
    println("Hoffman with $(length(w)) entries:")
    for k in sorted_keys
        c_str = coeff_to_str(w.terms[k])
        println("  ", lpad(c_str, pad), "  ", word_to_str(k))
    end
end

function sortedprint(idx::Index)
    if iszero(idx)
        println("Index (0)")
        return
    end

    coeffs = values(idx.terms)
    max_len = isempty(coeffs) ? 0 : maximum(c -> length(coeff_to_str(c)), coeffs)
    pad = max_len + 2
    
    # Sort by (weight, word)
    sorted_keys = sort(collect(keys(idx.terms)); by = t -> (sum(t.t), t.t))
    
    println("Index with $(length(idx)) entries:")
    for k in sorted_keys
        c_str = coeff_to_str(idx.terms[k])
        println("  ", lpad(c_str, pad), "  ", word_to_str(k))
    end
end

function sortedprint(r::Poly{A}) where {A <: Union{Hoffman, Index}}
    if iszero(r)
        println("Poly{$A} (0)")
        return
    end

    # Calculate padding
    pad = 0
    for ch in values(r.terms)
        if !isempty(ch.terms)
            coeffs = values(ch.terms)
            max_len = maximum(t -> length(coeff_to_str(t)), coeffs)
            pad = max(pad, max_len)
        end
    end
    pad += 2

    println("Poly{$A} with $(length(r.terms)) degree-entries:")
    
    degs = sort(collect(keys(r.terms)), rev=true)
    for d in degs
        ch = r.terms[d]
        if iszero(ch) continue end
        
        println("T^$d:")
        
        sorted_words = if A === Hoffman
            sort(collect(keys(ch.terms)); by = t -> (length(t), t.t))
        else # Index
            sort(collect(keys(ch.terms)); by = t -> (sum(t.t), t.t))
        end

        for wo in sorted_words
            c_str = coeff_to_str(ch.terms[wo])
            println("  ", lpad(c_str, pad), "  ", word_to_str(wo))
        end
    end
end

function sortedprint(r::Poly{Rational{BigInt}})
    if iszero(r)
        println("Poly{Rational{BigInt}} (0)")
        return
    end

    # Filter out zero-coefficient terms before calculating padding
    non_zero_terms = filter(p -> !iszero(p.second), r.terms)
    if isempty(non_zero_terms)
        println("Poly{Rational{BigInt}} (all coefficients are zero)")
        return
    end

    println("Poly{Rational{BigInt}} with $(length(non_zero_terms)) entries:")

    degs = sort(collect(keys(non_zero_terms)), rev=true)

    # Calculate padding for a neat table
    max_deg_len = maximum(d -> length(string(d)), keys(non_zero_terms))
    max_coeff_len = maximum(c -> length(coeff_to_str(c)), values(non_zero_terms))

    deg_pad = max(max_deg_len, length("Degree"))
    coeff_pad = max(max_coeff_len, length("Coefficient"))

    # Header
    println() # Add a blank line for spacing
    println(rpad("Degree", deg_pad), " | ", "Coefficient")
    println(repeat("─", deg_pad), "─┼─", repeat("─", coeff_pad))

    # Body
    for d in degs
        coeff = r.terms[d]
        println(rpad(string(d), deg_pad), " | ", coeff_to_str(coeff))
    end
end

"""
###################################################################################################
                                        GetProperty
###################################################################################################
"""

function getproperty(obj::Union{Hoffman, Index}, sym::Symbol)
    if sym === :sortshow
        sortedprint(obj)
        return nothing
    else
        return getfield(obj, sym)
    end
end

function getproperty(r::Poly, sym::Symbol)
    if sym === :sortshow
        sortedprint(r)
        return nothing
    elseif sym === :deg
        return isempty(r.terms) ? -1 : maximum(keys(r.terms))
    else
        return getfield(r, sym)
    end
end

function getproperty(w::Union{HoffmanWord, IndexWord}, sym::Symbol)
    if sym === :tovec
        return collect(Int, w.t)
    else
        return getfield(w, sym)
    end
end