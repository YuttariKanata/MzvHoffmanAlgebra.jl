using Test
using Random
using MzvHoffmanAlgebra

include("test_datas.jl")

# ==========================================
# Random Generators (from test_heavy.jl)
# ==========================================

function rand_hoffman_word(len_range)
    len = rand(len_range)
    return HoffmanWord(rand(0:1, len))
end

function rand_index_word(len_range, val_range)
    len = rand(len_range)
    return IndexWord(rand(val_range, len))
end

function rand_hoffman(num_terms, word_len_range)
    h = zero(Hoffman)
    for _ in 1:num_terms
        w = rand_hoffman_word(word_len_range)
        c = rand(-100:100) // rand(1:100)
        if !iszero(c)
            h = h + Hoffman(w, c)
        end
    end
    return h
end

function rand_index(num_terms, word_len_range, val_range)
    i = zero(Index)
    for _ in 1:num_terms
        w = rand_index_word(word_len_range, val_range)
        c = rand(-100:100) // rand(1:100)
        if !iszero(c)
            i = i + Index(w, c)
        end
    end
    return i
end

function rand_poly(deg_range, coeff_gen)
    deg = rand(deg_range)
    # Start with a zero poly. We need to know the type, so we generate one coeff to get the zero.
    dummy = coeff_gen()
    p = Poly(zero(dummy))
    
    for d in 0:deg
        if rand(Bool)
            c = coeff_gen()
            # Create term c * T^d
            # We use shift_degree on a constant poly to make T^d term
            term = shift_degree(Poly(c), d)
            p = p + term
        end
    end
    return p
end

@testset "MzvHoffmanAlgebra All Tests" begin

    # ==========================================
    # Comprehensive Unit Tests
    # ==========================================
    @testset "Unit Tests" begin
        
        @testset "1. Operators and Normalization" begin
            # Constructors
            op_up = OpUp(2)
            op_down = OpDown(3)
            @test op_up.cnt == 2
            @test op_down.cnt == 3
            
            # String conversion
            @test Operator("↑", 2) == Operator([OpUp(2)])
            @test Operator("↓", 1) == Operator([OpDown(1)])
            @test Operator("∂", 1) == Operator([OpDeriv(1,1)])
            @test Operator("∂2", 1) == Operator([OpDeriv(1,2)])
            
            # Normalization (clean)
            # OpDown -> OpUp(-cnt)
            op_clean = Operator([OpDown(2)])
            @test op_clean.ops[1] isa OpUp
            @test op_clean.ops[1].cnt == -2
            
            # Merging
            op_merge = Operator([OpUp(1), OpUp(2)])
            @test length(op_merge.ops) == 1
            @test op_merge.ops[1].cnt == 3
            
            # Zero count removal
            op_zero = Operator([OpUp(1), OpUp(-1)])
            @test isempty(op_zero.ops)
            @test isone(op_zero)
            
            # Multiplication (Concatenation)
            op1 = Operator("↑")
            op2 = Operator("←")
            op3 = op1 * op2
            @test length(op3.ops) == 2
            
            # Equality
            @test Operator("↑") == Operator("↑")
            @test Operator("↑") != Operator("↓")
        end

        @testset "2. Base Functions" begin
            # multinomial
            @test multinomial([1, 2]) == 3 # 3! / (1! 2!) = 6/2 = 3
            @test multinomial([2, 2]) == 6 # 4! / (2! 2!) = 24/4 = 6
            
            # idxprs / idxdprs
            # (1,1,0,0,1,0) -> [1, 3, 2]
            w = HoffmanWord(1, 1, 0, 0, 1, 0)
            @test idxprs(w) == [1, 3, 2]
            @test idxdprs([1, 3, 2]) == w
            
            # idxprs_r / idxdprs_r
            # (0,1,0,0,1,1) -> [2, 3, 1]
            w_r = HoffmanWord(0, 1, 0, 0, 1, 1)
            @test idxprs_r(w_r) == [2, 3, 1]
            @test idxdprs_r([2, 3, 1]) == w_r
        end

        @testset "3. Words and Predicates" begin
            # HoffmanWord
            w = HoffmanWord(1, 0)
            @test is_hoffmanword(w)
            @test !is_indexword(w)
            @test is_monomial(w)
            @test !iszero(w)
            @test !isone(w)
            @test isone(HoffmanWord())
            
            # IndexWord
            i = IndexWord(2, 1)
            @test is_indexword(i)
            @test !is_hoffmanword(i)
            
            # x, y
            @test x == HoffmanWord(0)
            @test y == HoffmanWord(1)
            
            # index()
            @test index(2, 1) == IndexWord(2, 1)
            
            # Admissibility
            # Left: ends with 0 (x) or >1
            # Hoffman Left: must start with 1 and end with 0. (y ... x)
            @test is_admissible(Hoffman(HoffmanWord(1, 0)), orientation=:left)
            @test !is_admissible(Hoffman(HoffmanWord(1)), orientation=:left)
            @test !is_admissible(Hoffman(HoffmanWord(0, 1)), orientation=:left)
            
            # Index Left: ends with > 1
            @test is_admissible(Index(IndexWord(2)), orientation=:left)
            @test !is_admissible(Index(IndexWord(1)), orientation=:left)
            
            # Right: starts with 0 (x) and ends with 1 (y)
            @test is_admissible(Hoffman(HoffmanWord(0, 1)), orientation=:right)
            @test !is_admissible(Hoffman(HoffmanWord(1, 0)), orientation=:right)
            
            # Index Right: starts with > 1
            @test is_admissible(Index(IndexWord(2)), orientation=:right)
            @test !is_admissible(Index(IndexWord(1)), orientation=:right)
        end

        @testset "4. Algebraic Elements (Hoffman/Index)" begin
            # Creation
            h = Hoffman(x)
            @test length(h) == 1
            @test h[x] == 1
            
            idx = Index(2)
            @test length(idx) == 1
            @test idx[index(2)] == 1
            
            # Arithmetic
            h2 = h + h
            @test h2[x] == 2
            h3 = h2 - h
            @test h3 == h
            
            h_prod = h * h # x * x = xx
            @test h_prod[x*x] == 1
            
            # Scalar
            @test (h * 2)[x] == 2
            @test (2 * h)[x] == 2
            @test (h // 2)[x] == 1//2
            
            # Power
            @test (h^3)[x*x*x] == 1
            
            # Zero/One
            @test iszero(zero(Hoffman))
            @test isone(one(Hoffman))
            @test zero(Hoffman) + h == h
            @test one(Hoffman) * h == h
        end

        @testset "5. Polynomials (Poly)" begin
            # Creation
            p = Poly(x)
            @test p[0] == Hoffman(x)
            
            # T
            @test T[1] == 1
            
            # Arithmetic
            p2 = p + T
            @test p2[0] == Hoffman(x)
            @test p2[1] == 1
            
            p3 = p * T
            @test p3[1] == Hoffman(x)
            
            # Shift degree
            p_shift = shift_degree(p, 2)
            @test p_shift[2] == Hoffman(x)
            
            # Power
            p_sq = (p + T)^2 # (x + T)^2 = x^2 + 2xT + T^2
            @test p_sq[0] == Hoffman(x*x)
            @test p_sq[1] == Hoffman(x) * 2
            @test p_sq[2] == 1
        end

        @testset "6. Conversions" begin
            # Hoffman <-> Index
            # Left: z_k = y x^{k-1} = 1 0^{k-1}
            # z_2 = 10
            h_z2 = Hoffman(HoffmanWord(1, 0))
            i_z2 = Index(h_z2, orientation=:left)
            @test i_z2[index(2)] == 1
            
            h_back = Hoffman(i_z2, orientation=:left)
            @test h_back == h_z2
            
            # Right: z_k = x^{k-1} y = 0^{k-1} 1
            # z_2 = 01
            h_z2_r = Hoffman(HoffmanWord(0, 1))
            i_z2_r = Index(h_z2_r, orientation=:right)
            @test i_z2_r[index(2)] == 1
            
            h_back_r = Hoffman(i_z2_r, orientation=:right)
            @test h_back_r == h_z2_r
            
            # Poly conversions
            p_h_valid = Poly(HoffmanWord(1, 0)) # z2
            p_i_valid = convert(Poly{Index}, p_h_valid, orientation=:left)
            @test p_i_valid[0] == Index(2)

            # New Hoffman constructors
            # Hoffman from Vector
            h_vec = Hoffman([0, 1]) # x y
            @test h_vec == Hoffman(x * y)
            
            # Hoffman from IndexWord with orientation
            # z_2 = yx (left) = 10
            iw = IndexWord(2)
            h_from_idx_left = Hoffman(iw, orientation=:left)
            @test h_from_idx_left == Hoffman(HoffmanWord(1, 0))
            
            # z_2 = xy (right) = 01
            h_from_idx_right = Hoffman(iw, orientation=:right)
            @test h_from_idx_right == Hoffman(HoffmanWord(0, 1))

            @testset "Poly Constructors" begin
                # Poly{Hoffman}
                p_h = Poly{Hoffman}(x)
                @test p_h isa Poly{Hoffman}
                @test p_h[0] == Hoffman(x)
                
                p_h2 = Poly{Hoffman}(Index(2), orientation=:left)
                @test p_h2[0] == Hoffman(HoffmanWord(1,0))

                # Poly{Index}
                p_i = Poly{Index}(index(2))
                @test p_i isa Poly{Index}
                @test p_i[0] == Index(2)
                
                p_i2 = Poly{Index}(Hoffman(x*y), orientation=:right) # xy = 01 -> 2
                @test p_i2[0] == Index(2)
            end
        end
        
        @testset "7. Edge Cases and Errors" begin
            # Empty
            @test HoffmanWord() * HoffmanWord() == HoffmanWord()
            
            # Invalid conversions
            @test_throws DomainError IndexWord(HoffmanWord(0), orientation=:left)
            
            # Negative power
            @test_throws ArgumentError Hoffman(x)^-1
            @test_throws DomainError Poly(x)^-1
        end
    end

    # ==========================================
    # Heavy Random Tests
    # ==========================================
    @testset "Heavy Random Tests" begin
        Random.seed!(42) # Ensure reproducibility

        @testset "1. Hoffman Algebra Laws (Stress Test)" begin
            for i in 1:100
                a = rand_hoffman(10, 0:8)
                b = rand_hoffman(10, 0:8)
                c = rand_hoffman(10, 0:8)

                # Associativity (+)
                @test (a + b) + c == a + (b + c)
                # Commutativity (+)
                @test a + b == b + a
                # Distributivity (*)
                @test a * (b + c) == a * b + a * c
                @test (a + b) * c == a * c + b * c
                # Associativity (*)
                @test (a * b) * c == a * (b * c)
            end
        end

        @testset "2. Index Algebra Laws (Stress Test)" begin
            for i in 1:100
                a = rand_index(10, 0:5, 1:5)
                b = rand_index(10, 0:5, 1:5)
                c = rand_index(10, 0:5, 1:5)

                @test (a + b) + c == a + (b + c)
                @test a + b == b + a
                @test a * (b + c) == a * b + a * c
                @test (a + b) * c == a * c + b * c
                @test (a * b) * c == a * (b * c)
            end
        end

        @testset "3. Conversion Roundtrip (Index -> Hoffman -> Index)" begin
            # Index -> Hoffman is always valid.
            # Hoffman -> Index is valid if the Hoffman element came from an Index element.
            for i in 1:50
                idx = rand_index(10, 1:6, 1:6)
                
                # Left Orientation
                h_left = Hoffman(idx, orientation=:left)
                idx_back_left = Index(h_left, orientation=:left)
                @test idx == idx_back_left

                # Right Orientation
                h_right = Hoffman(idx, orientation=:right)
                idx_back_right = Index(h_right, orientation=:right)
                @test idx == idx_back_right
            end
        end

        @testset "4. Poly Arithmetic & Mixed Types" begin
            for i in 1:50
                # Poly{Hoffman} laws
                gen_h = () -> rand_hoffman(5, 0:4)
                p1 = rand_poly(0:5, gen_h)
                p2 = rand_poly(0:5, gen_h)
                p3 = rand_poly(0:5, gen_h)

                @test (p1 + p2) + p3 == p1 + (p2 + p3)
                @test p1 * (p2 + p3) == p1 * p2 + p1 * p3

                # Mixed: Poly{Hoffman} + Poly{Rational}
                # T is Poly{Rational}
                h_const = rand_hoffman(1, 0:2)
                p_mixed = p1 + T # Should promote to Poly{Hoffman}
                @test typeof(p_mixed) == Poly{Hoffman}
                
                # Check specific value logic
                p_simple = Poly(h_const) + T
                @test p_simple[0] == h_const
                @test p_simple[1] == 1
            end
        end

        @testset "5. Regression & Edge Cases" begin
            # StackOverflow check for x + y
            @test (x + y) isa Hoffman
            @test length(x + y) == 2
            
            # Zero and One
            @test zero(Hoffman) + x == Hoffman(x)
            @test one(Hoffman) * x == Hoffman(x)
            @test x * one(Hoffman) == Hoffman(x)
            
            # Large Power
            w = HoffmanWord(1, 0)
            @test length(w^10) == 20
        end
    end

    @testset "8. Products (Shuffle & Stuffle)" begin
        # Shuffle Product
        # x ш y = xy + yx
        @test shuffle_product(Hoffman(x), Hoffman(y)) == Hoffman(x*y + y*x)
        # x ш x = 2xx
        @test shuffle_product(Hoffman(x), Hoffman(x)) == Hoffman(x*x * 2)
        # y ш y = 2yy
        @test shuffle_product(Hoffman(y), Hoffman(y)) == Hoffman(y*y * 2)

        # Stuffle Product
        # (1) * (1) = (1,1) + (1,1) + (2) = 2(1,1) + (2)
        i1 = Index(1)
        st = stuffle_product(i1, i1)
        expected = Index(1,1)*2 + Index(2)
        @test st == expected

        # Star Stuffle Product
        # (1) ⋆ (1) = (1,1) + (1,1) - (2) = 2(1,1) - (2)
        st_star = star_stuffle_product(i1, i1)
        expected_star = Index(1,1)*2 - Index(2)
        @test st_star == expected_star
    end

    @testset "9. Products with Index Type and Powers" begin
        # Cross-type products
        # Index shuffle: (1) ш (1) -> y ш y = 2yy -> 2(1,1)
        i1 = Index(1)
        sh_idx = shuffle_product(i1, i1, orientation=:left)
        @test sh_idx == Index(1,1)*2

        # Hoffman stuffle: y * y -> (1) * (1) = 2(1,1) + (2) -> 2yy + yx (left: z2=yx)
        h_y = Hoffman(y)
        st_hof = stuffle_product(h_y, h_y, orientation=:left)
        # (1)*(1) = 2(1,1) + (2). 
        # Left: (1,1)->yy, (2)->yx
        @test st_hof == Hoffman(y*y*2 + y*x)

        # Powers
        # y ш^2 = 2yy
        @test shuffle_pow(h_y, 2) == Hoffman(y*y*2)
        # (1) *^2 = 2(1,1) + (2)
        @test stuffle_pow(i1, 2) == Index(1,1)*2 + Index(2)
    end

    @testset "Optimized Index Shuffle" begin
        # Check if orientation is accepted (even if ignored by logic, it shouldn't error)
        i1 = Index(1)
        i2 = Index(2)
        # (1) ш (2) = 2*(1,2) + (2,1)
        expected = Index( shuffle_product(Hoffman(i1), Hoffman(i2)) )
        @test shuffle_product(i1, i2, orientation=:left) == expected
        @test shuffle_product(i1, i2, orientation=:right) == expected
    end

    @testset "10. Power Functions (power_by_squaring)" begin
        # Hoffman power
        h = x + y
        @test h^0 == one(Hoffman)
        @test h^1 == h
        @test h^2 == h*h
        @test h^3 == h*h*h
        @test h^4 == (h*h)*(h*h)

        # Index power
        i = index(2) + index(3)
        @test i^0 == one(Index)
        @test i^2 == i*i
        @test i^5 == i*i*i*i*i

        # Poly power
        p = Poly(x) + T
        @test p^0 == one(Poly{Hoffman})
        @test p^2 == p*p
        @test p^3 == p*p*p
        @test p^13 == foldl(*, fill(p, 13))
    end

    @testset "11. Special Power Functions" begin
        # st_index1_pow(2) = (1)*(1) = 2(1,1) + (2)
        # 2! / (1!1!) = 2. 2! / 2! = 1.
        st2 = st_index1_pow(2)
        @test st2[index(1,1)] == 2
        @test st2[index(2)] == 1
        
        # sh_y_pow(2) = y ш y = 2yy
        sh2 = sh_y_pow(2)
        @test sh2[HoffmanWord(1,1)] == 2
        
        # sh_y_pow(3) = 6yyy
        sh3 = sh_y_pow(3)
        @test sh3[HoffmanWord(1,1,1)] == 6
    end

    @testset "12. Regularization Map (rho)" begin
        # rho(T^2) = T^2 + zeta(2)
        # rho_t(2) should return (2, [])=>1 and (0, [2])=>1
        rt2 = rho_t(2)
        @test rt2[(2, Int[])] == 1
        @test rt2[(0, [2])] == 1
        
        # rho(T^4) = T^4 + 6zeta(2)T^2 - 8zeta(3)T + 6zeta(4) + 3zeta(2)^2
        rt4 = rho_t(4)
        @test rt4[(4, Int[])] == 1
        @test rt4[(2, [2])] == 6
        @test rt4[(1, [3])] == -8
        @test rt4[(0, [4])] == 6
        @test rt4[(0, [2,2])] == 3
    end

    @testset "13. Advanced Consistency & Robustness" begin
        Random.seed!(1234)
        
        # Access internal functions for testing
        mod = MzvHoffmanAlgebra
        
        gen_hoffman(n) = rand_hoffman(n, 0:5)
        gen_index(n) = rand_index(n, 1:3, 1:4) # Avoid empty words for strict conversion tests

        @testset "Squaring Optimization Consistency" begin
            println("\nTesting squaring optimizations...")
            for cnt in 1:10
                print("\r$cnt")

                h = gen_hoffman(5)
                i = gen_index(5)
                
                # Valid Hoffman for Index conversion (must start with y for :left)
                h_valid = Hoffman(i, orientation=:left)

                # Hoffman Shuffle: square(h) vs h * h
                @test mod.shuffle_square(h) == shuffle_product(h, h)
                
                # Index Stuffle: square(i) vs i * i
                @test mod.stuffle_square(i) == stuffle_product(i, i)
                
                # Index Star Stuffle: square(i) vs i * i
                @test mod.star_stuffle_square(i) == star_stuffle_product(i, i)
                
                # Cross-type wrappers
                # Index Shuffle
                @test mod.shuffle_square(i, orientation=:left) == shuffle_product(i, i, orientation=:left)
                
                # Hoffman Stuffle
                @test mod.stuffle_square(h_valid, orientation=:left) == stuffle_product(h_valid, h_valid, orientation=:left)
            end
        end

        @testset "Power Function Consistency" begin
            println("\nTesting power function consistency...")
            # Test small powers to cover _iter_pow and _bin_pow logic boundaries
            for n in 0:3
                print("\r$n")
                h = gen_hoffman(3)
                i = gen_index(3)
                
                # Naive calculation (foldl) vs Optimized Power
                h_sh_n = iszero(n) ? one(Hoffman) : foldl(shuffle_product, fill(h, n))
                i_st_n = iszero(n) ? one(Index) : foldl(stuffle_product, fill(i, n))
                
                @test shuffle_pow(h, n) == h_sh_n
                @test stuffle_pow(i, n) == i_st_n
            end
            println("\n")
            # Test special powers vs Generic powers
            for n in 0:5
                print("\r$n")
                # sh_y_pow(n) == y^{shuffle n}
                y_h = Hoffman(y)
                @test sh_y_pow(n) == shuffle_pow(y_h, n)
                
                # st_index1_pow(n) == (1)^{stuffle n}
                one_i = Index(index(1))
                @test st_index1_pow(n) == stuffle_pow(one_i, n)
            end
        end

        @testset "Regularization Consistency" begin
            # Convergent words should have trivial regularization (poly = const = self)
            
            # Convergent Index (left): ends with > 1 (e.g., (2,3))
            conv_idx = Index(index(2, 3))
            reg_poly = stuffle_regularization_polynomial(conv_idx, orientation=:left)
            @test length(reg_poly.terms) == 1
            @test reg_poly[0] == conv_idx
            
            # Convergent Hoffman (left): ends with 0 (x) (e.g., yx = 10)
            conv_hof = Hoffman(HoffmanWord(1, 0)) 
            reg_poly_h = shuffle_regularization_polynomial(conv_hof, orientation=:left)
            @test length(reg_poly_h.terms) == 1
            @test reg_poly_h[0] == conv_hof
        end
        
        @testset "Rho Map Linearity" begin
            # rho(A + B) = rho(A) + rho(B)
            pi1 = Poly(index(2)) + T 
            pi2 = Poly(index(3)) + T*2
            
            r1 = rho(pi1)
            r2 = rho(pi2)
            r_sum = rho(pi1 + pi2)
            
            @test r_sum == r1 + r2
        end
    end

    @testset "14. Data-Driven Verification" begin
        @testset "Shuffle Product Data" begin
            for ((a, b), expected) in SHUFFLE_PRODUCT_TEST_DATA
                if expected !== nothing
                    @test shuffle_product(a, b) == expected
                end
            end
        end
        
        @testset "Stuffle Product Data" begin
            for ((a, b), expected) in STUFFLE_PRODUCT_TEST_DATA
                if expected !== nothing
                    @test stuffle_product(a, b) == expected
                end
            end
        end

        @testset "Shuffle Regularization Data (Left)" begin
            for (input, expected) in SHUFFLE_REGULARIZATION_TEST_DATA_LEFT
                if expected !== nothing
                    @test shuffle_regularization_polynomial(input, orientation=:left) == expected
                end
            end
        end

        @testset "Stuffle Regularization Data (Left)" begin
            for (input, expected) in STUFFLE_REGULARIZATION_TEST_DATA_LEFT
                if expected !== nothing
                    @test stuffle_regularization_polynomial(input, orientation=:left) == expected
                end
            end
        end

        @testset "Rho Map Data" begin
            for (n, expected) in RHO_TEST_DATA
                if expected !== nothing
                    # rho takes Poly{Index}, so we convert T^n (which is Poly{Rational} or Poly{Hoffman} depending on context)
                    @test rho(convert(Poly{Index}, T^n)) == expected
                end
            end
        end

        @testset "Duals and Homomorphisms" begin
            # Dual: x(0) <-> y(1) and reverse
            # dual(xy) = dual(01) = 01 = xy
            @test dual(Hoffman(x*y)) == Hoffman(x*y)
            # dual(x) = y
            @test dual(Hoffman(x)) == Hoffman(y)
            
            # Hoffman Dual
            # Left: yxy (101) -> y * swap(01) = y * 10 = yyx (110)
            w = HoffmanWord(1, 0, 1)
            @test Hoffman_dual(Hoffman(w), orientation=:left) == Hoffman(HoffmanWord(1, 1, 0))
            
            # Landen Dual
            # Left: yx (10) -> (x+y) * (-y) = -xy -yy
            # Landen_dual(yx)
            h_yx = Hoffman(y*x)
            expected_landen = -(Hoffman(x*y) + Hoffman(y*y))
            @test Landen_dual(h_yx, orientation=:left) == expected_landen
            
            # Starword to word
            # Left: yx (10) -> y * x = yx
            # Left: yy (11) -> y * (x+y) = yx + yy
            @test starword_to_word(Hoffman(y*x), orientation=:left) == Hoffman(y*x)
            @test starword_to_word(Hoffman(y*y), orientation=:left) == Hoffman(y*x + y*y)
            
            # Homomorphism
            # x -> 2x, y -> 3y
            img = [Hoffman(x)*2, Hoffman(y)*3]
            # hom(xy) = hom(x)hom(y) = 2x * 3y = 6xy
            @test Hoffman_hom(HoffmanWord(0, 1), img) == Hoffman(x*y) * 6
            
            # Optimized Index Dual Consistency
            gen_index(n) = rand_index(n, 1:3, 1:4)
            for _ in 1:10
                
                idx = gen_index(5)
                is_admissible(idx; orientation=:left) || continue

                # Dual
                d_idx = dual(idx, orientation=:left)
                d_hof = dual(Hoffman(idx, orientation=:left))
                @test Hoffman(d_idx, orientation=:left) == d_hof
                
                # Hoffman Dual
                hd_idx = Hoffman_dual(idx, orientation=:left)
                hd_hof = Hoffman_dual(Hoffman(idx, orientation=:left), orientation=:left)
                @test Hoffman(hd_idx, orientation=:left) == hd_hof
            end

            @testset "Derivations" begin
                # Test dell(x, 2) for left orientation, based on oldhoffman.jl logic
                # ∂_2(x) = y * (x+y) * x = yxx + yyx
                expected_d2x = Hoffman(y*x*x + y*y*x)
                @test dell(Hoffman(x), 2, orientation=:left) == expected_d2x
            end
        end
    end

    @testset "15. Operators Test" begin

        @testset "Basic Operators on IndexWord" begin
            w = IndexWord(2, 3) # (2, 3)

            # OpMinus (minus)
            # Left action: remove from head
            @test minus(1) * w == IndexWord(3)
            # Right action: remove from tail
            @test w * minus(1) == IndexWord(2)
            
            # Boundary case for OpMinus
            # Note: Current implementation returns empty word (1) if length is insufficient for IndexWord,
            # but returns 0 for Index.
            @test minus(2) * w == IndexWord() 
            @test minus(3) * w == IndexWord()

            # OpUp (up)
            # Left: increase first index
            @test up(1) * w == IndexWord(3, 3)
            @test up(2) * w == IndexWord(4, 3)
            # Right: increase last index
            @test w * up(1) == IndexWord(2, 4)

            # OpDown (down)
            # Left: decrease first index
            @test down(1) * w == IndexWord(1, 3)
            # Right: decrease last index
            @test w * down(1) == IndexWord(2, 2)

            # OpLeft (left) / OpRight (right)
            # Left: prepend 1s
            @test left(1) * w == IndexWord(1, 2, 3)
            @test left(2) * w == IndexWord(1, 1, 2, 3)
            # Right: append 1s
            @test w * right(1) == IndexWord(2, 3, 1)
            @test w * right(2) == IndexWord(2, 3, 1, 1)
        end

        @testset "Operators on Index" begin
            # Index = (2,3) + 2*(4)
            idx = Index(2,3) + Index(4)*2

            # OpMinus
            # (2,3) -> (3), (4) -> length 1 <= 1 -> removed (becomes 0)
            @test minus(1) * idx == Index(3)
            
            # OpUp
            # (2,3) -> (3,3), (4) -> (5)
            @test up(1) * idx == Index(3,3) + Index(5)*2

            # Linearity check
            idx2 = Index(3,1,4)
            @test up(1) * (idx + idx2) == (up(1) * idx) + (up(1) * idx2)
        end

        @testset "Special Operators" begin
            # OpTau (τ) - Dual
            # dual((3)) = (1, 2)
            @test τ(1) * IndexWord(3) == IndexWord(1, 2)
            # dual((1, 2)) = (3)
            @test τ(1) * IndexWord(1, 2) == IndexWord(3)
            # τ^2 = id
            @test τ(2) * IndexWord(2, 3) == IndexWord(2, 3)

            # OpDeriv (∂)
            # ∂_1 (1) = ∂_1(y). 
            # In left orientation: ∂_1(y) = -xyy - yxy = -(2,1) - (1,2) ? 
            # Let's check implementation logic via test.
            # dell(y, 1) -> xyn1 = y(x+y)^0 x = yx = (1,1). image=[yx, -yx].
            # d(y) = -yx = -(2).
            @test ∂(1) * IndexWord(1) == -Index(2)
            
            # d(x) = yx = (2). x is not directly representable in IndexWord unless we use Hoffman.
            # But IndexWord(2) = yx. d(yx) = d(y)x + yd(x) = (yx)x + y(-yx) = yyx - yxx
            # = (1,2) - (3)
            @test ∂(1) * IndexWord(2) == Index(1,2) - Index(3)
        end

        @testset "Operator Algebra" begin
            # Composition
            # op = Left(1) then Up(1)
            # w = (1). Left(1)->(1,1). Up(1)->(2,1).
            op = up(1) * left(1)
            @test op * IndexWord(1) == IndexWord(2, 1)
            
            # Associativity
            op2 = (up(1) * left(1)) * right(1)
            op3 = up(1) * (left(1) * right(1))
            @test op2.ops == op3.ops
            
            # Power
            # Up(1)^2 = Up(2)
            @test (up(1))^2 * IndexWord(1) == IndexWord(3)
            
            # Operator construction from Word
            # WordtoOperator((2,1))
            # Should be: Right(1) * Up(1) * Right(1) * Up(0)
            # = Right(1) * Up(1) * Right(1)
            op_w = WordtoOperator(IndexWord(2, 1))
            # Apply to empty word (1)
            # 1 * Right(1) -> (1)
            # (1) * Up(1) -> (2)
            # (2) * Right(1) -> (2,1)
            # (2,1) * Up(0) -> (2,1)
            # So WordtoOperator(w) * 1 should generate w
            @test IndexWord()* op_w == IndexWord(2, 1)
        end

        @testset "String Representation" begin
            op = up(1) * ∂(1)
            # Just check it runs without error
            io = IOBuffer()
            show(io, "text/plain", op)
            s = String(take!(io))
            @test length(s) > 0
        end
    end
end