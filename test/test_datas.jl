# ==========================================
# Test Data Collection
# ==========================================
# Please fill in the `nothing` values with the expected output.
# Example: (Hoffman(x), Hoffman(y)) => Hoffman(x*y + y*x)

# ------------------------------------------
# 1. Shuffle Product
# ------------------------------------------
const SHUFFLE_PRODUCT_TEST_DATA = Dict{Tuple{Hoffman, Hoffman}, Any}(
    # Weight 2
    (Hoffman(x), Hoffman(y)) => (x*y + y*x),
    
    # Weight 3
    (Hoffman(x*y), Hoffman(x)) => (x*y*x + 2*x*x*y),
    (Hoffman(y*y), Hoffman(y)) => (3*y*y*y),

    # Weight 4
    (Hoffman(x*y*x), Hoffman(y)) => (y*x*y*x + x*y*x*y + 2*x*y*y*x),
    (Hoffman(y*y*y), Hoffman(y)) => (4*y*y*y*y),
    
    # Weight 5 (Complex case)
    (Hoffman(x*y*x*y), Hoffman(y)) => (2*x*y*y*x*y + y*x*y*x*y + 2*x*y*x*y*y),

    (Hoffman(y*y*x+x*y), Hoffman(x*x+x*y*x*y)) => (x*y*x*y*y*x*y + x*y*x*x + y*x*x*y*x + 2*x*y*y*x*x 
                                                 + 3*x*y*y*y*x*y*x + 2*x*x*y*x + 4*x*y*y*x*y*y*x 
                                                 + 6*x*y*y*y*x*x*y + 3*x*y*x*y*x*y + 3*y*y*x*x*x 
                                                 + y*y*x*y*x*y*x + 3*x*y*y*x*y*x*y + 4*x*y*x*x*y*y 
                                                 + 4*x*x*y*y*x*y + 2*y*x*y*x*x + 2*y*x*y*y*x*y*x 
                                                 + 2*y*y*x*y*x*x*y + 4*x*x*y*x*y*y + 2*y*x*y*x*y*y*x 
                                                 + 4*y*x*y*y*x*x*y + 2*y*y*x*x*y*x*y + x*y*x*y*x 
                                                 + 2*y*x*y*x*y*x*y + 3*x*x*x*y + x*x*y*y*x + 3*x*y*x*y*y*y*x),
)

# ------------------------------------------
# 2. Stuffle Product
# ------------------------------------------
const STUFFLE_PRODUCT_TEST_DATA = Dict{Tuple{Index, Index}, Any}(
    # Weight 3
    (Index(1), Index(2)) => Index(1,2) + Index(2,1) + Index(3),
    (Index(1,1), Index(1)) => 3*Index(1,1,1) + Index(2,1) + Index(1,2),

    # Weight 4
    (Index(2,1), Index(1)) => Index(1,2,1) + Index(3,1) + Index(2,2) + 2*Index(2,1,1),
    (Index(1,1,1), Index(1)) => Index(1,2,1) + 4*Index(1,1,1,1) + Index(1,1,2) + Index(2,1,1),
    
    # Weight 5
    (Index(2,1,1), Index(1)) => Index(3,1,1) + Index(1,2,1,1) + Index(2,2,1) + Index(2,1,2) + 3*Index(2,1,1,1),

    (Index(3,1,4), Index(1,5,9)) => Index(3, 1, 6, 4, 9) + Index(1, 3, 5, 9, 1, 4) + 2Index(3, 1, 1, 5, 9, 4) 
                                  + Index(1, 3, 1, 4, 5, 9) + Index(3, 1, 5, 1, 9, 4) + Index(1, 5, 3, 10, 4) 
                                  + Index(1, 8, 1, 13) + Index(1, 5, 3, 1, 4, 9) + Index(4, 1, 5, 13) 
                                  + Index(1, 5, 3, 1, 9, 4) + Index(3, 1, 4, 1, 5, 9) + Index(1, 3, 6, 13) 
                                  + Index(3, 1, 5, 9, 1, 4) + 2Index(3, 1, 1, 4, 5, 9) + Index(1, 3, 1, 5, 13) 
                                  + Index(1, 3, 5, 1, 13) + Index(3, 2, 4, 5, 9) + Index(4, 5, 1, 9, 4) 
                                  + Index(4, 6, 4, 9) + Index(4, 5, 1, 4, 9) + Index(1, 5, 3, 9, 1, 4) 
                                  + Index(3, 2, 5, 9, 4) + Index(1, 5, 12, 1, 4) + Index(3, 2, 5, 4, 9) 
                                  + Index(3, 2, 9, 9) + Index(3, 1, 6, 13) + 2Index(3, 1, 1, 5, 13) + Index(4, 6, 13) 
                                  + Index(1, 5, 9, 3, 1, 4) + Index(3, 1, 5, 5, 9) + Index(3, 1, 5, 1, 13) 
                                  + Index(1, 8, 1, 9, 4) + Index(4, 5, 9, 1, 4) + Index(1, 8, 1, 4, 9) 
                                  + Index(4, 1, 4, 5, 9) + Index(1, 3, 1, 9, 9) + Index(1, 5, 3, 1, 13) 
                                  + Index(4, 1, 5, 9, 4) + Index(4, 6, 9, 4) + Index(1, 3, 5, 10, 4) 
                                  + Index(1, 3, 1, 5, 4, 9) + Index(4, 5, 10, 4) + Index(4, 1, 5, 4, 9) 
                                  + Index(1, 3, 5, 1, 4, 9) + Index(4, 1, 9, 9) + Index(1, 3, 6, 9, 4) 
                                  + Index(1, 3, 6, 4, 9) + Index(1, 3, 1, 5, 9, 4) + Index(1, 3, 5, 1, 9, 4) 
                                  + Index(1, 8, 9, 1, 4) + Index(4, 5, 1, 13) + 2Index(3, 1, 1, 9, 9) 
                                  + Index(3, 2, 5, 13) + Index(3, 1, 5, 10, 4) + 2Index(3, 1, 1, 5, 4, 9) 
                                  + Index(1, 8, 10, 4) + Index(3, 1, 5, 1, 4, 9) + Index(3, 1, 6, 9, 4),
)

# ------------------------------------------
# 3. Shuffle Regularization
# ------------------------------------------
const SHUFFLE_REGULARIZATION_TEST_DATA_LEFT = Dict{Hoffman, Any}(
    # y (T)
    Hoffman(y) => Poly{Hoffman}(T),
    # y^2 (T^2/2)
    Hoffman(y*y) => Poly{Hoffman}((1//2)*T^2),
    # yxy
    Hoffman(y*x*y) => Poly{Hoffman}(y*x*T - 2*y*y*x),
    # y^3
    Hoffman(y*y*y) => Poly{Hoffman}((1//6)*T^3),

    Hoffman(x*x*x*x*y*y*x*x*x*y*y*x*x*y*y*y) => 1//6*x*x*x*x*y*y*x*x*x*y*y*x*x*T^3 + ( - 1//2*x*x*x*x*y*y*x*y*x*x*y*y*x*x - 1//2*x*y*x*x*x*y*y*x*x*x*y*y*x*x - 1//2*x*x*x*y*x*y*y*x*x*x*y*y*x*x - 3//2*x*x*x*x*y*y*x*x*x*y*y*y*x*x - 1//2*x*x*y*x*x*y*y*x*x*x*y*y*x*x - 1//2*x*x*x*x*y*y*x*x*y*x*y*y*x*x - 1//2*x*x*x*x*y*y*x*x*x*y*y*x*y*x - 3//2*x*x*x*x*y*y*y*x*x*x*y*y*x*x - 1//2*y*x*x*x*x*y*y*x*x*x*y*y*x*x)T^2 + ( + x*y*x*x*x*y*y*x*x*y*x*y*y*x*x + x*x*x*y*x*y*y*x*y*x*x*y*y*x*x + 3*x*x*x*x*y*y*y*x*y*x*x*y*y*x*x + 3*x*y*x*x*x*y*y*y*x*x*x*y*y*x*x + y*x*y*x*x*x*y*y*x*x*x*y*y*x*x + x*x*x*y*x*y*y*x*x*x*y*y*x*y*x + x*y*x*y*x*x*y*y*x*x*x*y*y*x*x + 3*x*x*x*x*y*y*y*x*x*x*y*y*x*y*x + x*y*x*x*y*x*y*y*x*x*x*y*y*x*x + x*x*x*x*y*y*x*x*y*y*x*y*y*x*x + y*x*x*x*x*y*y*x*y*x*x*y*y*x*x + x*x*x*x*y*y*x*y*y*x*x*y*y*x*x + 3*x*x*x*y*x*y*y*x*x*x*y*y*y*x*x + 9*x*x*x*x*y*y*y*x*x*x*y*y*y*x*x + 3*x*x*x*x*y*y*x*x*x*y*y*y*x*y*x + y*x*x*x*x*y*y*x*x*x*y*y*x*y*x + y*y*x*x*x*x*y*y*x*x*x*y*y*x*x + x*x*x*x*y*y*x*x*x*y*y*x*y*y*x + x*x*x*x*y*y*x*y*x*x*y*y*x*y*x + 6*x*x*x*x*y*y*x*x*x*y*y*y*y*x*x + 3*y*x*x*x*x*y*y*x*x*x*y*y*y*x*x + x*x*y*x*x*y*y*x*y*x*x*y*y*x*x + 3*x*x*x*x*y*y*x*y*x*x*y*y*y*x*x + x*x*x*y*x*y*y*x*x*y*x*y*y*x*x + 3*x*x*x*x*y*y*y*x*x*y*x*y*y*x*x + x*x*y*x*x*y*y*x*x*x*y*y*x*y*x + x*y*y*x*x*x*y*y*x*x*x*y*y*x*x + 3*x*x*x*y*x*y*y*y*x*x*x*y*y*x*x + 6*x*x*x*x*y*y*y*y*x*x*x*y*y*x*x + 3*x*x*y*x*x*y*y*x*x*x*y*y*y*x*x + y*x*x*x*x*y*y*x*x*y*x*y*y*x*x + x*x*x*y*y*x*y*y*x*x*x*y*y*x*x + x*y*x*x*x*y*y*x*y*x*x*y*y*x*x + x*x*x*x*y*y*x*y*x*y*x*y*y*x*x + 3*y*x*x*x*x*y*y*y*x*x*x*y*y*x*x + x*y*x*x*x*y*y*x*x*x*y*y*x*y*x + y*x*x*y*x*x*y*y*x*x*x*y*y*x*x + y*x*x*x*y*x*y*y*x*x*x*y*y*x*x + x*x*y*x*x*y*y*x*x*y*x*y*y*x*x + 3*x*y*x*x*x*y*y*x*x*x*y*y*y*x*x + x*x*x*x*y*y*x*x*y*x*y*y*x*y*x + 3*x*x*y*x*x*y*y*y*x*x*x*y*y*x*x + x*x*y*y*x*x*y*y*x*x*x*y*y*x*x + 3*x*x*x*x*y*y*x*x*y*x*y*y*y*x*x + x*x*y*x*y*x*y*y*x*x*x*y*y*x*x)T - x*x*y*x*x*y*y*x*x*y*y*x*y*y*x*x - 3*x*x*y*x*x*y*y*y*x*x*x*y*y*x*y*x - 3*x*y*x*y*x*x*y*y*x*x*x*y*y*y*x*x - x*x*x*x*y*y*x*x*y*x*y*y*x*y*y*x - x*y*x*x*x*y*y*x*y*y*x*x*y*y*x*x - 3*x*x*x*x*y*y*x*x*y*y*x*y*y*y*x*x - y*y*x*x*x*y*x*y*y*x*x*x*y*y*x*x - x*x*x*y*y*x*y*y*x*y*x*x*y*y*x*x - 3*x*y*x*x*y*x*y*y*x*x*x*y*y*y*x*x - 3*y*x*x*x*x*y*y*x*y*x*x*y*y*y*x*x - 3*x*x*y*x*x*y*y*y*x*x*y*x*y*y*x*x - x*x*y*y*x*x*y*y*x*x*x*y*y*x*y*x - 3*y*x*x*x*x*y*y*y*x*y*x*x*y*y*x*x - y*y*x*y*x*x*x*y*y*x*x*x*y*y*x*x - x*x*y*x*y*x*y*y*x*x*x*y*y*x*y*x - 6*x*y*x*x*x*y*y*y*y*x*x*x*y*y*x*x - 6*x*x*x*y*x*y*y*x*x*x*y*y*y*y*x*x - 3*x*x*x*x*y*y*x*x*y*x*y*y*y*x*y*x - 3*y*x*y*x*x*x*y*y*y*x*x*x*y*y*x*x - 3*x*x*x*x*y*y*x*y*y*x*x*y*y*y*x*x - y*y*x*x*y*x*x*y*y*x*x*x*y*y*x*x - 18*x*x*x*x*y*y*y*x*x*x*y*y*y*y*x*x - x*x*y*y*x*x*y*y*x*x*y*x*y*y*x*x - 3*y*y*x*x*x*x*y*y*x*x*x*y*y*y*x*x - x*y*x*x*x*y*y*x*x*y*x*y*y*x*y*x - x*x*y*x*y*x*y*y*x*x*y*x*y*y*x*x - 3*x*y*x*y*x*x*y*y*y*x*x*x*y*y*x*x - y*x*x*y*x*x*y*y*x*y*x*x*y*y*x*x - x*x*x*y*x*y*y*x*y*x*x*y*y*x*y*x - y*x*x*x*y*x*y*y*x*y*x*x*y*y*x*x - 3*x*y*x*x*y*x*y*y*y*x*x*x*y*y*x*x - 10*x*x*x*x*y*y*x*x*x*y*y*y*y*y*x*x - y*x*y*x*x*x*y*y*x*x*x*y*y*x*y*x - 3*x*x*x*x*y*y*y*x*y*x*x*y*y*x*y*x - x*y*x*x*x*y*y*x*x*y*y*x*y*y*x*x - x*y*y*x*x*y*x*y*y*x*x*x*y*y*x*x - 3*x*y*x*x*x*y*y*y*x*x*x*y*y*x*y*x - 6*y*x*x*x*x*y*y*x*x*x*y*y*y*y*x*x - x*x*x*y*x*y*y*x*y*x*y*x*y*y*x*x - 3*x*x*y*x*x*y*y*x*y*x*x*y*y*y*x*x - 6*x*x*x*x*y*y*x*y*x*x*y*y*y*y*x*x - x*y*y*y*x*x*x*y*y*x*x*x*y*y*x*x - y*x*y*x*x*x*y*y*x*x*y*x*y*y*x*x - 3*x*x*x*x*y*y*y*x*y*x*y*x*y*y*x*x - 3*x*y*x*x*x*y*y*y*x*x*y*x*y*y*x*x - x*y*x*y*x*x*y*y*x*x*x*y*y*x*y*x - x*x*x*y*x*y*y*x*x*x*y*y*x*y*y*x - 3*x*x*y*x*x*y*y*y*x*y*x*x*y*y*x*x - x*y*y*x*y*x*x*y*y*x*x*x*y*y*x*x - 3*y*y*x*x*x*x*y*y*y*x*x*x*y*y*x*x - x*y*x*x*y*x*y*y*x*x*x*y*y*x*y*x - x*x*x*x*y*y*x*x*y*y*x*y*y*x*y*x - 3*x*x*x*y*x*y*y*x*x*y*x*y*y*y*x*x - 3*x*x*x*x*y*y*y*x*x*x*y*y*x*y*y*x - x*x*x*y*y*y*x*y*y*x*x*x*y*y*x*x - 3*x*y*y*x*x*x*y*y*x*x*x*y*y*y*x*x - x*x*y*y*x*x*y*y*x*y*x*x*y*y*x*x - y*x*x*x*x*y*y*x*y*x*x*y*y*x*y*x - 9*x*x*x*x*y*y*y*x*x*y*x*y*y*y*x*x - x*y*x*y*x*x*y*y*x*x*y*x*y*y*x*x - 3*x*x*x*y*x*y*y*x*x*x*y*y*y*x*y*x - x*x*x*x*y*y*x*x*y*y*y*x*y*y*x*x - x*y*x*x*y*x*y*y*x*x*y*x*y*y*x*x - x*x*x*x*y*y*x*y*y*x*x*y*y*x*y*x - x*x*y*x*y*x*y*y*x*y*x*x*y*y*x*x - 9*x*x*x*y*x*y*y*y*x*x*x*y*y*y*x*x - 9*x*x*x*x*y*y*y*x*x*x*y*y*y*x*y*x - 18*x*x*x*x*y*y*y*y*x*x*x*y*y*y*x*x - y*x*x*x*x*y*y*x*y*x*y*x*y*y*x*x - 3*x*x*x*x*y*y*x*x*x*y*y*y*x*y*y*x - y*y*y*x*x*x*x*y*y*x*x*x*y*y*x*x - x*x*x*x*y*y*x*x*x*y*y*x*y*y*y*x - 6*x*x*y*x*x*y*y*x*x*x*y*y*y*y*x*x - y*y*x*x*x*x*y*y*x*x*x*y*y*x*y*x - x*x*x*x*y*y*x*y*y*x*y*x*y*y*x*x - y*x*x*x*x*y*y*x*x*x*y*y*x*y*y*x - 3*x*x*x*y*y*x*y*y*x*x*x*y*y*y*x*x - y*x*x*y*x*y*x*y*y*x*x*x*y*y*x*x - 3*y*x*x*x*x*y*y*x*x*y*x*y*y*y*x*x - x*x*x*y*x*y*y*x*y*y*x*x*y*y*x*x - x*x*x*x*y*y*x*y*x*x*y*y*x*y*y*x - 3*x*y*x*x*x*y*y*x*y*x*x*y*y*y*x*x - 3*x*y*y*x*x*x*y*y*y*x*x*x*y*y*x*x - y*y*x*x*x*x*y*y*x*x*y*x*y*y*x*x - 6*x*x*x*x*y*y*x*x*x*y*y*y*y*x*y*x - y*x*x*x*y*y*x*y*y*x*x*x*y*y*x*x - 3*x*x*x*x*y*y*x*y*x*y*x*y*y*y*x*x - 3*x*x*x*x*y*y*y*x*y*y*x*x*y*y*x*x - y*x*y*x*x*x*y*y*x*y*x*x*y*y*x*x - 3*x*y*x*x*x*y*y*y*x*y*x*x*y*y*x*x - 3*y*x*x*x*x*y*y*x*x*x*y*y*y*x*y*x - 6*x*x*x*y*x*y*y*y*y*x*x*x*y*y*x*x - 9*y*x*x*x*x*y*y*y*x*x*x*y*y*y*x*x - x*x*y*x*x*y*y*x*y*x*x*y*y*x*y*x - 10*x*x*x*x*y*y*y*y*y*x*x*x*y*y*x*x - 3*x*x*x*x*y*y*x*y*x*x*y*y*y*x*y*x - y*x*x*y*y*x*x*y*y*x*x*x*y*y*x*x - x*y*x*y*x*x*y*y*x*y*x*x*y*y*x*x - x*x*x*y*x*y*y*x*x*y*x*y*y*x*y*x - 3*y*x*x*y*x*x*y*y*x*x*x*y*y*y*x*x - x*x*y*x*x*y*y*x*y*x*y*x*y*y*x*x - x*y*x*x*y*x*y*y*x*y*x*x*y*y*x*x - 3*x*x*x*y*y*x*y*y*y*x*x*x*y*y*x*x - 3*x*x*x*x*y*y*y*x*x*y*x*y*y*x*y*x - x*y*y*x*x*x*y*y*x*x*x*y*y*x*y*x - y*x*x*x*x*y*y*x*y*y*x*x*y*y*x*x - 3*y*x*x*x*y*x*y*y*x*x*x*y*y*y*x*x - x*x*y*x*x*y*y*x*x*x*y*y*x*y*y*x - 6*x*y*x*x*x*y*y*x*x*x*y*y*y*y*x*x - x*x*x*x*y*y*x*y*y*y*x*x*y*y*x*x - x*x*x*y*x*y*y*x*x*y*y*x*y*y*x*x - x*x*y*y*x*y*x*y*y*x*x*x*y*y*x*x - 3*x*x*y*x*x*y*y*x*x*y*x*y*y*y*x*x - 3*x*x*x*y*x*y*y*y*x*x*x*y*y*x*y*x - 6*y*x*x*x*x*y*y*y*y*x*x*x*y*y*x*x - 6*x*x*x*x*y*y*y*y*x*x*x*y*y*x*y*x - x*x*y*x*y*y*x*y*y*x*x*x*y*y*x*x - 3*x*x*x*x*y*y*y*x*x*y*y*x*y*y*x*x - x*y*y*x*x*x*y*y*x*x*y*x*y*y*x*x - 3*x*x*y*x*x*y*y*x*x*x*y*y*y*x*y*x - y*y*x*x*x*x*y*y*x*y*x*x*y*y*x*x - 9*x*x*y*x*x*y*y*y*x*x*x*y*y*y*x*x - 3*x*x*x*y*x*y*y*y*x*x*y*x*y*y*x*x - x*x*y*y*y*x*x*y*y*x*x*x*y*y*x*x - x*x*x*y*y*x*y*y*x*x*x*y*y*x*y*x - 6*x*x*x*x*y*y*y*y*x*x*y*x*y*y*x*x - y*x*x*x*x*y*y*x*x*y*x*y*y*x*y*x - 3*y*x*x*y*x*x*y*y*y*x*x*x*y*y*x*x - x*y*x*x*x*y*y*x*y*x*x*y*y*x*y*x - x*x*x*x*y*y*x*y*x*y*x*y*y*x*y*x - 3*x*x*y*y*x*x*y*y*x*x*x*y*y*y*x*x - 3*y*x*x*x*y*x*y*y*y*x*x*x*y*y*x*x - 6*x*x*x*x*y*y*x*x*y*x*y*y*y*y*x*x - x*x*x*y*y*x*y*y*x*x*y*x*y*y*x*x - x*x*y*x*x*y*y*x*y*y*x*x*y*y*x*x - y*x*y*x*x*y*x*y*y*x*x*x*y*y*x*x - 3*y*x*x*x*x*y*y*y*x*x*x*y*y*x*y*x - y*x*x*x*x*y*y*x*x*y*y*x*y*y*x*x - 3*x*x*y*x*y*x*y*y*x*x*x*y*y*y*x*x - x*y*x*x*x*y*y*x*y*x*y*x*y*y*x*x - x*x*x*x*y*y*x*y*x*y*y*x*y*y*x*x - x*y*x*x*x*y*y*x*x*x*y*y*x*y*y*x - 6*x*x*y*x*x*y*y*y*y*x*x*x*y*y*x*x - 3*y*x*x*x*x*y*y*y*x*x*y*x*y*y*x*x - y*x*x*y*x*x*y*y*x*x*x*y*y*x*y*x - y*x*y*y*x*x*x*y*y*x*x*x*y*y*x*x - 3*x*y*x*x*x*y*y*x*x*y*x*y*y*y*x*x - x*y*x*y*x*y*x*y*y*x*x*x*y*y*x*x - x*y*y*x*x*x*y*y*x*y*x*x*y*y*x*x - x*y*x*x*y*y*x*y*y*x*x*x*y*y*x*x - y*x*x*x*y*x*y*y*x*x*x*y*y*x*y*x - y*x*y*x*y*x*x*y*y*x*x*x*y*y*x*x - 3*x*x*x*y*x*y*y*x*y*x*x*y*y*y*x*x - 3*x*y*x*x*x*y*y*x*x*x*y*y*y*x*y*x - y*x*x*y*x*x*y*y*x*x*y*x*y*y*x*x - x*x*y*x*x*y*y*x*x*y*x*y*y*x*y*x - 9*x*x*x*x*y*y*y*x*y*x*x*y*y*y*x*x - 3*y*x*y*x*x*x*y*y*x*x*x*y*y*y*x*x - 3*x*x*y*y*x*x*y*y*y*x*x*x*y*y*x*x - 9*x*y*x*x*x*y*y*y*x*x*x*y*y*y*x*x - 3*x*x*x*y*x*y*y*y*x*y*x*x*y*y*x*x - x*y*x*y*y*x*x*y*y*x*x*x*y*y*x*x - 6*x*x*x*x*y*y*y*y*x*y*x*x*y*y*x*x - y*x*x*x*y*x*y*y*x*x*y*x*y*y*x*x - 3*x*x*y*x*y*x*y*y*y*x*x*x*y*y*x*x,
)

# ------------------------------------------
# 4. Stuffle Regularization
# ------------------------------------------
const STUFFLE_REGULARIZATION_TEST_DATA_LEFT = Dict{Index, Any}(
    # (1) (T)
    Index(1) => Poly{Index}(T),
    # (1,1)
    Index(1,1) => Poly{Index}(1//2*T^2 - 1//2*Index(2)),
    # (1,2)
    Index(1,2) => Poly{Index}(Index(1,2)),
    # (1,1,1)
    Index(1,1,1) => Poly{Index}(1//6*T^3 - 1//2*Index(2)*T + 1//3*Index(3)),

    6*Index(1,4,1,3,1,1,1) => Index(1, 4, 1, 3)T^3 
                       + ( - 6*Index(1, 1, 4, 1, 3) - 3*Index(2, 4, 1, 3) - 3*Index(1, 4, 2, 3) - 6*Index(1, 4, 1, 1, 3) - 3*Index(1, 5, 1, 3) - 3*Index(1, 4, 1, 4))T^2 
                       + ( + 3*Index(1, 4, 3, 3) + 6*Index(2, 4, 1, 4) + 9*Index(2, 1, 4, 1, 3) + 3*Index(3, 4, 1, 3) + 3*Index(1, 4, 1, 5) + 12*Index(1, 5, 1, 1, 3) 
                       + 12*Index(1, 4, 1, 1, 4) + 9*Index(1, 4, 2, 1, 3) + 3*Index(1, 6, 1, 3) + 24*Index(1, 1, 4, 1, 1, 3) + 12*Index(1, 1, 4, 2, 3) 
                       + 12*Index(2, 4, 1, 1, 3) + 6*Index(1, 5, 2, 3) + 12*Index(1, 1, 5, 1, 3) - 3*Index(1, 4, 1, 3, 2) + 18*Index(1, 1, 1, 4, 1, 3) 
                       + 9*Index(1, 2, 4, 1, 3) + 6*Index(1, 5, 1, 4) + 6*Index(2, 5, 1, 3) + 12*Index(1, 1, 4, 1, 4) + 6*Index(1, 4, 2, 4) + 6*Index(2, 4, 2, 3) 
                       + 18*Index(1, 4, 1, 1, 1, 3) + 9*Index(1, 4, 1, 2, 3))T 
                       - 24*Index(1, 1, 1, 1, 4, 1, 3) + 3*Index(1, 4, 2, 3, 2) - 6*Index(3, 4, 1, 1, 3) - 6*Index(1, 1, 4, 3, 3) - 9*Index(1, 5, 1, 2, 3) 
                       - 12*Index(1, 4, 2, 1, 1, 3) - 6*Index(1, 4, 1, 1, 5) - 9*Index(1, 4, 1, 2, 4) - 6*Index(2, 5, 1, 4) - 18*Index(1, 1, 1, 5, 1, 3) 
                       - Index(4, 4, 1, 3) - 24*Index(1, 1, 4, 1, 1, 4) - 9*Index(2, 4, 1, 2, 3) - 24*Index(1, 1, 5, 1, 1, 3) - Index(1, 4, 1, 6) - Index(1, 7, 1, 3) 
                       - 12*Index(1, 5, 1, 1, 4) - 36*Index(1, 1, 4, 1, 1, 1, 3) - 9*Index(1, 2, 4, 2, 3) - 18*Index(1, 1, 1, 4, 1, 4) - Index(1, 4, 4, 3) 
                       - 12*Index(2, 4, 1, 1, 4) - 12*Index(2, 5, 1, 1, 3) - 3*Index(1, 6, 1, 4) - 36*Index(1, 1, 1, 4, 1, 1, 3) - 9*Index(1, 2, 4, 1, 4) 
                       + 6*Index(1, 1, 4, 1, 3, 2) - 6*Index(1, 5, 2, 4) - 6*Index(1, 4, 2, 2, 3) + 3*Index(1, 4, 1, 4, 2) - 6*Index(2, 2, 4, 1, 3) - 18*Index(1, 1, 4, 2, 1, 3) 
                       - 6*Index(1, 6, 1, 1, 3) - 6*Index(2, 4, 2, 4) - 3*Index(3, 5, 1, 3) - 18*Index(1, 1, 4, 1, 2, 3) - 9*Index(1, 5, 2, 1, 3) - 9*Index(1, 4, 2, 1, 4) 
                       - 2*Index(1, 4, 1, 3, 3) - 3*Index(1, 4, 3, 4) - 12*Index(1, 1, 5, 2, 3) - 4*Index(1, 3, 4, 1, 3) - 18*Index(1, 1, 1, 4, 2, 3) - 9*Index(2, 4, 2, 1, 3) 
                       - 18*Index(1, 4, 1, 1, 1, 4) - 12*Index(1, 1, 5, 1, 4) - 18*Index(2, 1, 4, 1, 1, 3) - 4*Index(1, 4, 3, 1, 3) - 24*Index(1, 4, 1, 1, 1, 1, 3) 
                       - 6*Index(1, 1, 6, 1, 3) - 9*Index(2, 1, 5, 1, 3) - 12*Index(2, 1, 1, 4, 1, 3) - 3*Index(1, 5, 1, 5) - 3*Index(3, 4, 2, 3) + 6*Index(1, 4, 1, 1, 3, 2) 
                       - 4*Index(3, 1, 4, 1, 3) - 3*Index(2, 4, 1, 5) - 12*Index(1, 1, 4, 2, 4) - 6*Index(1, 1, 4, 1, 5) - 12*Index(1, 4, 1, 2, 1, 3) - 12*Index(1, 4, 1, 1, 2, 3) 
                       - 3*Index(2, 6, 1, 3) + 3*Index(1, 5, 1, 3, 2) - 12*Index(1, 1, 2, 4, 1, 3) - 9*Index(2, 1, 4, 2, 3) - 18*Index(1, 5, 1, 1, 1, 3) 
                       - 6*Index(2, 5, 2, 3) - 3*Index(1, 5, 3, 3) + 3*Index(2, 4, 1, 3, 2) - 9*Index(2, 1, 4, 1, 4) - 18*Index(2, 4, 1, 1, 1, 3) - 3*Index(1, 4, 2, 5) 
                       - 3*Index(2, 4, 3, 3) - 18*Index(1, 2, 4, 1, 1, 3) - 3*Index(1, 6, 2, 3) - 3*Index(3, 4, 1, 4) - 9*Index(1, 2, 5, 1, 3) - 12*Index(1, 2, 1, 4, 1, 3)
)

# ------------------------------------------
# 5. Rho Map (T^n -> Poly)
# ------------------------------------------
# Input is integer n (degree of T)
const RHO_TEST_DATA = Dict{Int, Any}(
    2 => T^2 + y*x,
    4 => T^4 + 6*y*x*T^2 - 8*y*x*x*T + 6*y*x*y*x + 12*y*y*x*x + 6*y*x*x*x,
    6 => T^6 + 15*y*x*T^4 - 40*y*x*x*T^3 + ( + 90*y*x*y*x + 180*y*y*x*x + 90*y*x*x*x)*T^2 + ( - 180*y*x*y*x*x - 60*y*x*x*y*x - 360*y*y*x*x*x - 144*y*x*x*x*x)*T + 170*y*x*x*y*x*x + 840*y*y*x*x*x*x + 120*y*x*x*x*x*x + 180*y*y*x*x*y*x + 180*y*x*y*y*x*x + 420*y*x*y*x*x*x + 45*y*x*x*x*y*x + 540*y*y*y*x*x*x + 90*y*x*y*x*y*x + 360*y*y*x*y*x*x,
)

# ------------------------------------------
# 6. Derivation (Dell)
# ------------------------------------------
# Input is (Hoffman, n)
const DELL_TEST_DATA = Dict{Tuple{Hoffman, Int}, Any}(
    (Hoffman(x), 1) => nothing,
    (Hoffman(x*y), 1) => nothing,
    (Hoffman(y), 1) => nothing,
    (Hoffman(x), 2) => nothing,
    (Hoffman(x*y*x), 1) => nothing,
)