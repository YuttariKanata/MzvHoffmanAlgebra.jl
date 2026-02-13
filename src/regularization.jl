#[ regularization.jl ]#

# This file defines functions for regularization and the rho map.

"""
###################################################################################################
                                        Regularization
###################################################################################################
"""

# Helper to split a word into (remainder, count) of consecutive values at the "divergent" end.
# If orientation is :left, divergence is at the end (suffix).
# If orientation is :right, divergence is at the start (prefix).
function _split_consecutive(w::W, val::Int, orientation::Symbol) where W <: AbstractWord
    t = w.t
    len = length(t)
    cnt = 0
    
    if orientation === :left
        # Count from end
        while cnt < len && t[end-cnt] == val
            cnt += 1
        end
        return (W(t[1:end-cnt]), cnt)
    elseif orientation === :right
        # Count from start
        while cnt < len && t[1+cnt] == val
            cnt += 1
        end
        return (W(t[cnt+1:end]), cnt)
    else
        throw(ArgumentError("orientation must be :left or :right"))
    end
end

# Generic regularization kernel
"""
    _regularize(input::P, product_op, pow_op, val::Int, orientation::Symbol) where P <: Union{Hoffman, Index} -> Poly{P}

 - product_op: 積の定義(`shuffle_product` または `stuffle_product`).
 - pow_op: `T^n` の生成ルール.
 - シャッフル正規化なら sh_y_pow (`y^{⧢ n} = n! y^n`).
 - 調和正規化なら st_index1_pow (`(1)^{*n}`).
 - val: 発散判定に使う値 (`Hoffman`なら `1`、`Index`なら `1`) .
 - orientation: 左右どちらを見るか.
"""
function _regularize(input::P, product_op, pow_op, val::Int, orientation::Symbol) where P <: Union{Hoffman, Index}
    
    z = copy(input)
    res = zero(Poly{P})
    
    while !iszero(z)
        # 1. Find max consecutive count
        max_cnt = -1
        for w in keys(z.terms)
            _, cnt = _split_consecutive(w, val, orientation)
            if cnt > max_cnt
                max_cnt = cnt
            end
        end
        
        # If max_cnt is 0, it means no divergence found (or we are done with divergent parts).
        # We still need to add these terms to T^0.
        
        # 2. Collect terms with max_cnt
        targets = Vector{Tuple{keytype(P), valtype(P)}}()
        for (w, c) in z.terms
            _, cnt = _split_consecutive(w, val, orientation)
            if cnt == max_cnt
                push!(targets, (w, c))
            end
        end
        
        # 3. Process targets
        pow_term = pow_op(max_cnt) # T^k equivalent in the algebra
        
        for (w, c) in targets
            rem, _ = _split_consecutive(w, val, orientation)
            
            # Construct the term: rem * val^k (or val^k * rem depending on orientation/product)
            # Since shuffle/stuffle are commutative, order doesn't strictly matter for the result,
            # but we must be consistent with how we split.
            # We split `w` into `rem` and `val^k`.
            # So `w` approx `rem * val^k`.
            correction = product_op(P(rem), pow_term)
            
            denom = correction[w]
            if iszero(denom)
                error("Algorithm failure: Leading term not found in product.")
            end
            
            diff = c // denom
            
            # Subtract from z
            z = z - correction * diff
            
            # Add to result polynomial
            current_poly_term = get(res.terms, max_cnt, zero(P))
            res.terms[max_cnt] = current_poly_term + P(rem) * diff
        end
    end
    
    return res
end

"""
    stuffle_regularization_polynomial(x; orientation::Symbol=:left)

Calculates the regularization polynomial of `x` with respect to the stuffle product.
This treats the divergent index `1` as a variable `T`.
Returns a `Poly{Index}` where the coefficient of `T^0` is the regularized value.
"""
function stuffle_regularization_polynomial(idx::Index; orientation::Symbol=:left)::Poly{Index}
    _regularize(idx, stuffle_product, st_index1_pow, 1, orientation)
end

"""
    stuffle_regularization_polynomial(h::Hoffman; orientation::Symbol=:left) -> Poly{Hoffman}

Calculates the regularization polynomial of `h` with respect to the stuffle product.
Converts to Index, computes the polynomial, and converts back to Poly{Hoffman}.
"""
function stuffle_regularization_polynomial(h::Hoffman; orientation::Symbol=:left)::Poly{Hoffman}
    idx = Index(h, orientation=orientation)
    poly_i = stuffle_regularization_polynomial(idx, orientation=orientation)
    return Poly{Hoffman}(poly_i, orientation=orientation)
end
stuffle_regularization_polynomial(w::IndexWord; orientation::Symbol=:left) = stuffle_regularization_polynomial(Index(w), orientation=orientation)
stuffle_regularization_polynomial(w::HoffmanWord; orientation::Symbol=:left) = stuffle_regularization_polynomial(Hoffman(w), orientation=orientation)

"""
    shuffle_regularization_polynomial(x; orientation::Symbol=:left)

Calculates the regularization polynomial of `x` with respect to the shuffle product.
This treats the divergent letter `y` (1) as a variable `T`.
Returns a `Poly{Hoffman}` where the coefficient of `T^0` is the regularized value.
"""
function shuffle_regularization_polynomial(h::Hoffman; orientation::Symbol=:left)::Poly{Hoffman}
    _regularize(h, shuffle_product, sh_y_pow, 1, orientation)
end

"""
    shuffle_regularization_polynomial(idx::Index; orientation::Symbol=:left) -> Poly{Index}

Calculates the regularization polynomial of `idx` with respect to the shuffle product.
Converts to Hoffman, computes the polynomial, and converts back to Poly{Index}.
"""
function shuffle_regularization_polynomial(idx::Index; orientation::Symbol=:left)::Poly{Index}
    h = Hoffman(idx, orientation=orientation)
    poly_h = shuffle_regularization_polynomial(h, orientation=orientation)
    return Poly{Index}(poly_h, orientation=orientation)
end
shuffle_regularization_polynomial(w::HoffmanWord; orientation::Symbol=:left) = shuffle_regularization_polynomial(Hoffman(w), orientation=orientation)
shuffle_regularization_polynomial(w::IndexWord; orientation::Symbol=:left) = shuffle_regularization_polynomial(Index(w), orientation=orientation)

"""
    reg_st(x; orientation::Symbol=:left)

Returns the constant term of the stuffle regularization polynomial of `x`.
This is the "finite part" or "regularized value" with respect to the stuffle product.
"""
reg_st(x::Union{Index, IndexWord}; kw...) = get(stuffle_regularization_polynomial(x; kw...).terms, 0, zero(Index))
reg_st(h::Union{Hoffman, HoffmanWord}; kw...) = Hoffman(reg_st(Index(h;kw... )); kw... )

"""
    reg_sh(x; orientation::Symbol=:left)

Returns the constant term of the shuffle regularization polynomial of `x`.
This is the "finite part" or "regularized value" with respect to the shuffle product.
"""
reg_sh(x::Union{Hoffman, HoffmanWord}; kw...) = get(shuffle_regularization_polynomial(x; kw...).terms, 0, zero(Hoffman))
reg_sh(x::Union{Index, IndexWord}; kw...) = Index(reg_sh(Hoffman(x; kw...)); kw...)

"""
###################################################################################################
                                        Regularization Map (rho)
###################################################################################################
"""

struct CombSum
    a::Vector{Int}
    s::Int
end

function sum_combinations(n::Int, k::Int)::Vector{CombSum}
    if k == 0
        return [CombSum(Int[], 0)]
    end

    results = Vector{CombSum}()
    a = fill(2, k)
    S = 2 * k

    if S > n
        return results
    end

    while true
        push!(results, CombSum(copy(a), S))

        advanced = false
        suffix_sum = 0

        for i in k:-1:1
            suffix_sum += a[i]
            ai = a[i] + 1
            rest = k - i
            
            new_sum = (S - suffix_sum) + ai * (rest + 1)

            if new_sum <= n
                S = new_sum
                a[i] = ai
                for j in i+1:k
                    a[j] = ai
                end
                advanced = true
                break
            end
        end

        if !advanced
            break
        end
    end

    return results
end

"""
    rho_t(m::Int)

Computes the map ρ(T^m) as a dictionary mapping (degree, zeta_indices) to coefficient.
ρ(e^{Tu}) = A(u)e^{Tu} where A(u) = exp(sum_{n>=2} (-1)^n/n * zeta(n) * u^n).
"""
function rho_t(m::Int)::Dict{Tuple{Int, Vector{Int}}, Rational{BigInt}}
    result = Dict{Tuple{Int, Vector{Int}}, Rational{BigInt}}()
    mfact = factorial(BigInt(m))

    for k in 0:(m >> 1)
        m_k = mfact // factorial(BigInt(k))
        for cs in sum_combinations(m, k)
            j = m - cs.s
            
            zetas = cs.a
            coeff = m_k // factorial(BigInt(j))
            
            for n in zetas
                coeff *= (-1)^n // n
            end
            
            key = (j, zetas)
            result[key] = get(result, key, 0//1) + coeff
        end
    end
    return result
end

"""
    rho(poly::Poly{Index}; orientation::Symbol=:left) -> Poly{Hoffman}

Applies the map `ρ` to a polynomial `poly`.
The map `ρ` is defined by `ρ(e^{Tu}) = A(u)e^{Tu}` where `A(u) = exp(sum_{n>=2} (-1)^n/n * zeta(n) * u^n)`.
It maps `T^m` to a linear combination of `T` and zeta values.
This map connects the stuffle regularization parameter `T` to the shuffle regularization parameter `T`.
"""
function rho(ri::Poly{Index}; orientation::Symbol=:left)::Poly{Hoffman}
    rh = zero(Poly{Hoffman})
    for (deg, coeff) in ri.terms
        coeff_h = Hoffman(coeff, orientation=orientation)
        rT = rho_t(deg)
        for ((rhodeg, zetas), ratcoeff) in rT
            term_val = coeff_h
            for n in zetas
                z_n = Hoffman(IndexWord(n), orientation=orientation)
                term_val = shuffle_product(term_val, z_n)
            end
            current = get(rh.terms, rhodeg, zero(Hoffman))
            rh.terms[rhodeg] = current + term_val * ratcoeff
        end
    end
    return rh
end