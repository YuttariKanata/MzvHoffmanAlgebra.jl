using MzvHoffmanAlgebra
using Test

@testset "MzvHoffmanAlgebra.jl" begin
    # Write your tests here.

    @testset "types.jl" begin
        ################################### types.jl ###################################
        @testset "Operators" begin
            op1 = OpMinus(123)
            @test typeof(op1) === OpMinus
            @test typeof(op1) <: AbstractOp
            @test op1.cnt == 123
            op1 = OpEta()
            @test typeof(op1) === OpEta
            @test op1.cnt == 1
            op2 = OpDeriv(19)
            @test typeof(op2) === OpDeriv
            @test op2.cnt == 19 && op2.n == 1
            op3 = OpDeriv(17,13)
            @test op3.cnt == 13 && op3.n == 17
            op4 = OpDeriv()
            @test op4.cnt == 1 && op4.n == 1
            @test typeof(Operator()) === Operator
        end
        
        @testset "Hoffman MZV" begin
            @test 0 isa NN
            @test -13//17 isa NN
            @test Word((1,2,3)).t == (1,2,3)
    
            @test isdefined(MzvHoffmanAlgebra,:_INDEX_ORIENTATION)
            @test_nowarn get_index_orientation()
            now_index_orientation = get_index_orientation()
            @test_nowarn set_index_orientation!(true)
            @test get_index_orientation() == true
            @test_nowarn set_index_orientation!(false)
            @test get_index_orientation() == false
            @test_nowarn set_index_orientation!(now_index_orientation)

            @test isdefined(MzvHoffmanAlgebra,:_OMIT_COUNTS)
            
            @test_nowarn a=Poly{Rational{BigInt}}()
            @test_nowarn a=Poly{Rational{BigInt}}(Dict{Int64,Rational{BigInt}}(1=>0,2=>-19//17))
            a = Poly{Rational{BigInt}}(Dict{Int64,Rational{BigInt}}(1=>0,2=>-19//17))
            @test a.terms[1] == 0 && a.terms[2] == -19//17

            @test isdefined(MzvHoffmanAlgebra,:T)
            @test length(T.terms) == 1 && T.terms[1] == 1
        end
    end

    
end
;