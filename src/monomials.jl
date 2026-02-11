#[ monomials.jl ]#

# This file defines low-level monomial product algorithms.

"""
###################################################################################################
                                        Monomial Products
###################################################################################################
"""

# Helper to add terms to a dictionary
function _add!(d::Dict{W, C}, w::W, c::C) where {W, C}
    iszero(c) && return
    new_val = get(d, w, zero(C)) + c
    if iszero(new_val)
        delete!(d, w)
    else
        d[w] = new_val
    end
end

# Implements the Suffix-recursive definition of shuffle product:
# 1 ш w = w ш 1 = w
# wu ш w'u' = (wu ш w')u' + (w ш w'u')u
# This is algebraically equivalent to the Prefix-recursive definition but more efficient for Vectors (appending).
function monomial_sh(a::Vector{Int},b::Vector{Int})::Vector{Vector{Int}}
    lb = lastindex(b)
    la = lastindex(a)
    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[1:i]]
    end

    for i in 1:lb
        # Previous row's h1[1] was [Int[]] (empty word).
        # New h1[1] should be [b[1:i]].
        # We need to update h1 in place or use a buffer. 
        # The logic below iterates j from 1 to la, updating h1[j+1].
        # We need to be careful not to overwrite h1[j] before using it for h1[j+1] calculation?
        # Actually, the DP relation is:
        # sh(a[1:j], b[1:i]) = sh(a[1:j], b[1:i-1]) * b[i] + sh(a[1:j-1], b[1:i]) * a[j]
        # Let H(j, i) = sh(a[1:j], b[1:i])
        # H(j, i) depends on H(j, i-1) [current h1[j+1]] and H(j-1, i) [newly computed h1[j]]
        
        # We can update in place if we store the "new h1[j]" in a temporary variable or iterate carefully.
        # But here, let's just use a full new row to be safe and clear, then swap.
        
        next_h1 = Vector{Vector{Vector{Int}}}(undef, la+1)
        next_h1[1] = [b[1:i]]

        for j in 1:la
            # term1: h1[j+1] (which is H(j, i-1)) appended with b[i]
            # term2: next_h1[j] (which is H(j-1, i)) appended with a[j]
            
            list1 = h1[j+1]
            list2 = next_h1[j]
            
            len1 = length(list1)
            len2 = length(list2)
            
            new_terms = Vector{Vector{Int}}(undef, len1 + len2)
            
            # Manual loop to avoid intermediate allocations from vcat broadcasting
            idx = 1
            for w in list1
                new_terms[idx] = push!(copy(w), b[i])
                idx += 1
            end
            for w in list2
                new_terms[idx] = push!(copy(w), a[j])
                idx += 1
            end
            
            next_h1[j+1] = new_terms
        end
        h1 = next_h1
    end

    return h1[end]
end

# Implements the Suffix-recursive definition of stuffle product.
function monomial_st(a::Vector{Int},b::Vector{Int})::Vector{Vector{Int}}
    lb = lastindex(b)
    la = lastindex(a)
    h1 = Vector{Vector{Vector{Int}}}(undef,la+1)

    h1[1] = [Int[]]
    for i in 1:la
        h1[i+1] = [a[1:i]]
    end

    for i in 1:lb
        # H(j, i) = H(j, i-1)*b[i] + H(j-1, i)*a[j] + H(j-1, i-1)*(a[j]+b[i])
        # h1[j+1] currently holds H(j, i-1)
        # We need H(j-1, i) which we are computing in this row.
        # We need H(j-1, i-1) which is h1[j] from the previous iteration (before update).
        
        # To avoid complex dependency tracking, let's use next_h1
        next_h1 = Vector{Vector{Vector{Int}}}(undef, la+1)
        next_h1[1] = [b[1:i]]
        
        for j in 1:la
            list_prev_i   = h1[j+1]      # H(j, i-1)
            list_prev_j   = next_h1[j]   # H(j-1, i)
            list_prev_all = h1[j]        # H(j-1, i-1)
            
            len1 = length(list_prev_i)
            len2 = length(list_prev_j)
            len3 = length(list_prev_all)
            
            new_terms = Vector{Vector{Int}}(undef, len1 + len2 + len3)
            idx = 1
            
            # term1: H(j, i-1) * b[i]
            for w in list_prev_i
                new_terms[idx] = push!(copy(w), b[i])
                idx += 1
            end
            
            # term2: H(j-1, i) * a[j]
            for w in list_prev_j
                new_terms[idx] = push!(copy(w), a[j])
                idx += 1
            end
            
            # term3: H(j-1, i-1) * (a[j]+b[i])
            sum_val = a[j] + b[i]
            for w in list_prev_all
                new_terms[idx] = push!(copy(w), sum_val)
                idx += 1
            end
            
            next_h1[j+1] = new_terms
        end
        h1 = next_h1
    end

    return h1[end]
end

# Type alias for the dictionary used in Index
const TermDict = Dict{IndexWord, Rational{BigInt}}
"""
    monomial_sh_i(aI::Vector{Int}, bI::Vector{Int}) -> Index

Computes the shuffle product of two index words (given as vectors of integers) using a block-based dynamic programming approach.
This method is significantly faster and more memory-efficient than the standard recursive definition for Index words.
"""
function monomial_sh_i(aI::Vector{Int}, bI::Vector{Int})::Index
    
    # 1. Initialize state vectors with prefixes
    # Instead of Vector{Index}, we use Vector{TermDict} to avoid immutable struct overhead
    total_len_a = sum(aI)
    ai = Vector{TermDict}(undef, total_len_a)
    
    idx = 1
    current_prefix = Int[]
    for len in aI
        for k in 1:len
            # Create prefix: (current_prefix..., k)
            # We can optimize this by pushing to a buffer, but for initialization it's fine
            term = [current_prefix; k]
            d = TermDict()
            d[IndexWord(term)] = one(Rational{BigInt})
            ai[idx] = d
            idx += 1
        end
        push!(current_prefix, len)
    end

    total_len_b = sum(bI)
    bi = Vector{TermDict}(undef, total_len_b)
    
    idx = 1
    current_prefix = Int[]
    for len in bI
        for k in 1:len
            term = [current_prefix; k]
            d = TermDict()
            d[IndexWord(term)] = one(Rational{BigInt})
            bi[idx] = d
            idx += 1
        end
        push!(current_prefix, len)
    end

    # 2. Block-based DP
    idxa = 0
    idxb = 0
    
    # Helper to add terms with digit appending
    function _append_and_add!(dest::TermDict, src::TermDict, digit::Int, scale::BigInt)
        iszero(scale) && return
        for (w, c) in src
            new_t = (w.t..., digit)
            new_w = IndexWord(new_t)
            # Optimized dictionary update
            new_val = get(dest, new_w, zero(Rational{BigInt})) + c * scale
            if !iszero(new_val)
                dest[new_w] = new_val
            else
                delete!(dest, new_w)
            end
        end
    end

    for len_a in aI
        # Process block A
        for len_b in bI
            # Process block B
            
            # We need to update the segment of `ai` corresponding to the current block of A
            # and the segment of `bi` corresponding to the current block of B.
            # The update depends on the *current* values of these segments.
            
            # To avoid creating new dictionaries for every cell, we can update in place or use buffers.
            # Since `ai` and `bi` segments are read and written, we need to be careful.
            # The logic in al.jl uses `vshv!` which writes to `ra` and `rb` based on `a` and `b`.
            # Here `a` and `b` are the current states of the segments.
            
            # Extract views of the current segments
            range_a = (idxa+1):(idxa+len_a)
            range_b = (idxb+1):(idxb+len_b)
            
            # Create new dictionaries for the next state to avoid mutation issues during read
            next_ai_segment = [TermDict() for _ in 1:len_a]
            next_bi_segment = [TermDict() for _ in 1:len_b]
            
            # Calculate next states
            # Update A segment
            for k in 1:len_a
                s1 = len_b + k - 1
                # From A
                for i in 1:k
                    scale = binomial(BigInt(s1-i), BigInt(len_b-1))
                    _append_and_add!(next_ai_segment[k], ai[range_a[i]], s1 - i + 1, scale)
                end
                # From B
                for i in 1:len_b
                    scale = binomial(BigInt(s1-i), BigInt(k-1))
                    _append_and_add!(next_ai_segment[k], bi[range_b[i]], s1 - i + 1, scale)
                end
            end
            
            # Update B segment
            for k in 1:len_b-1
                s1 = len_a + k - 1
                # From A
                for i in 1:len_a
                    scale = binomial(BigInt(s1-i), BigInt(k-1))
                    _append_and_add!(next_bi_segment[k], ai[range_a[i]], s1 - i + 1, scale)
                end
                # From B
                for i in 1:k
                    scale = binomial(BigInt(s1-i), BigInt(len_a-1))
                    _append_and_add!(next_bi_segment[k], bi[range_b[i]], s1 - i + 1, scale)
                end
            end
            # Optimization: The last element of B's segment is the same as the last element of A's segment
            next_bi_segment[len_b] = copy(next_ai_segment[len_a]) 

            # Apply updates
            for k in 1:len_a
                ai[range_a[k]] = next_ai_segment[k]
            end
            for k in 1:len_b
                bi[range_b[k]] = next_bi_segment[k]
            end
            
            idxb += len_b
        end
        idxa += len_a
        idxb = 0
    end

    return Index(ai[end])
end