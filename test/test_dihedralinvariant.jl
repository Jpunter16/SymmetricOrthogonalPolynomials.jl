using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials, MultivariateOrthogonalPolynomials, BlockArrays, StaticArrays, Test
using SymmetricOrthogonalPolynomials: dihedralconversion

@testset "DihedralInvariant" begin
    Q = DihedralInvariant(Legendre())
    P² = KronPolynomial(Legendre(), Legendre())
    x,y = 0.1,0.2
    for K = Block.(1:10)
        @test Q[SVector(x,y),K] == Q[SVector(y,x),K] == Q[SVector(x,-y),K]
    end


    f = (x,y) -> exp(x^2+y^2)*cos(x^2 * y^2)

    c = transform(P², splat(f))
    N = (Int(last(BlockArrays.blockcolsupport(c)))+1)÷2
    R = dihedralconversion(N)
    d = R * c[Block.(1:2:2N)]
   
    
    @test Q[SVector(x,y),Block.(1:N)]'d ≈ f(x,y)

    @testset "mass matrix" begin
        N= 4
        R = dihedralconversion(N)
        @test R*grammatrix(P²)[Block.(1:2:2N),Block.(1:2:2N)]*R' ≈ grammatrix(Q)[Block.(1:N), Block.(1:N)] ≈ (Q'Q)[Block.(1:N), Block.(1:N)]
    end

    @testset "variable coefficient" begin
        X = jacobimatrix(Val(1), P²)
        Y = jacobimatrix(Val(2), P²)
        N = 20
        R = dihedralconversion(N)
        d = R * c[Block.(1:2:2N)]
        A = R * (X^2 + Y^2)[Block.(1:2:2N), Block.(1:2:2N)] * R'
        B = R * (X^2 * Y^2)[Block.(1:2:2N), Block.(1:2:2N)] * R'
        C = R*((X^2 - Y^2)^2)[Block.(1:2:2N), Block.(1:2:2N)] * R'

        @test Q[SVector(x,y),axes(A,1)]'A*d ≈ (x^2 + y^2) * f(x,y)
        @test Q[SVector(x,y),axes(B,1)]'B*d ≈ (x^2 * y^2) * f(x,y)
        @test Q[SVector(x,y),axes(C,1)]'C*d ≈ (x^2 - y^2)^2 * f(x,y)

        Ã = grammatrix(Q)[Block.(1:N),Block.(1:N)] * A
        B̃ = grammatrix(Q)[Block.(1:N),Block.(1:N)] * B
        C̃ = grammatrix(Q)[Block.(1:N),Block.(1:N)] * C
        @test Ã ≈ Ã'
        @test B̃ ≈ B̃'
        @test C̃ ≈ C̃'
    end

    @testset "laplacian" begin
        W = Weighted(Ultraspherical(3/2))
        Q = DihedralInvariant(W)
        W² = KronPolynomial(W, W)
        N = 10
        R = dihedralconversion(N)
        M = R * grammatrix(W²)[Block.(1:2:2N), Block.(1:2:2N)] * R'
        M₁ = grammatrix(W)

        @test M₁[1,1]^2 ≈ M[1,1]
        @test 2M₁[1,1]M₁[1,3]/sqrt(2) ≈ M[1,2]
        @test 0 ≈ M[1,3]
        @test M₁[1,3]^2 ≈ M[1,4]

        # (2k,2j) = (2,0)
        # (2k̃,2j̃) = (2,0)
        @test M₁[3,3]M₁[1,1] + M₁[1,3]^2 ≈ M[2,2]
        # (2k,2j) = (2,0)
        # (2k̃,2j̃) = (4,0)
        @test M₁[1,1]M₁[3,5] ≈ M[Block(2,3)[1,1]]
        # (2k,2j) = (2,0)
        # (2k̃,2j̃) = (2,2)
        @test 2M₁[3,3]M₁[1,3]/sqrt(2) ≈ M[Block(2,3)[1,2]]
        # (2k,2j) = (2,0)
        # (2k̃,2j̃) = (4,2)
        @test M₁[3,5]M₁[1,3] ≈ M[Block(2,4)[1,2]]

        # (2k,2j) = (4,0)
        # (2k̃,2j̃) = (4,0)
        @test M₁[5,5]M₁[1,1] ≈ M[Block(3,3)[1,1]]
        # (2k,2j) = (4,0)
        # (2k̃,2j̃) = (2,2)
        @test 2M₁[5,3]M₁[1,3]/sqrt(2) ≈ M[Block(3,3)[1,2]]
        # (2k,2j) = (4,0)
        # (2k̃,2j̃) = (6,0)
        @test M₁[5,7]M₁[1,1] ≈ M[Block(3,4)[1,1]]
        # (2k,2j) = (4,0)
        # (2k̃,2j̃) = (4,2)
        @test M₁[5,5]M₁[1,3] ≈ M[Block(3,4)[1,2]]
        # (2k,2j) = (4,0)
        # (2k̃,2j̃) = (6,2)
        @test M₁[5,7]M₁[1,3] ≈ M[Block(3,5)[1,2]]


        # (2k,2j) = (2,2)
        # (2k̃,2j̃) = (2,2)
        @test M₁[3,3]^2 ≈ M[Block(3,3)[2,2]]
        # (2k,2j) = (2,2)
        # (2k̃,2j̃) = (6,0)
        @test M₁[3,7]M₁[3,1] ≈ M[Block(3,4)[2,1]] ≈ 0
        # (2k,2j) = (2,2)
        # (2k̃,2j̃) = (4,2)
        @test 2/sqrt(2) * M₁[3,5]M₁[3,3] ≈ M[Block(3,4)[2,2]]
        # (2k,2j) = (2,2)
        # (2k̃,2j̃) = (4,4)
        @test M₁[3,5]^2 ≈ M[Block(3,5)[2,3]]


        for K = 1:2:9, J = 1:2:9
            for k=1:(K+1)÷2 , j = 1:(J+1)÷2
                if (j == (J+1)÷2 && k == (K+1)÷2) || (j ≠ (J+1)÷2 && k ≠ (K+1)÷2)
                    @test M₁[2K-2k+1,2J-2j+1]M₁[2k-1,2j-1] ≈ M[Block(K,J)[k,j]]
                else
                    @test 2/sqrt(2) * M₁[2K-2k+1,2J-2j+1]M₁[2k-1,2j-1] ≈ M[Block(K,J)[k,j]]
                end
            end
        end

        for K = 1:2:9, J = 2:2:9
            for k=1:(K+1)÷2 , j = 1:(J+1)÷2
                if k == (K+1)÷2
                    @test 2/sqrt(2) * M₁[2K-2k+1,2J-2j+1]M₁[2k-1,2j-1] ≈ M[Block(K,J)[k,j]]
                else
                    @test M₁[2K-2k+1,2J-2j+1]M₁[2k-1,2j-1] ≈ M[Block(K,J)[k,j]]
                end
            end
        end

        for K = 2:2:11, J = 2:2:11
            for k=1:(K+1)÷2, j = 1:(J+1)÷2
                if (j == (J+1)÷2 && k == (K+1)÷2)
                    @test M₁[2K-2k+1,2J-2j+1]M₁[2k-1,2j-1] + M₁[2K-2k+1,2j-1]M₁[2J-2j+1,2k-1] ≈ M[Block(K,J)[k,j]]
                else
                    @test M₁[2K-2k+1,2J-2j+1]M₁[2k-1,2j-1] ≈ M[Block(K,J)[k,j]]
                end
            end
        end

        for K = 2:2:9, J = 1:2:9
            for k=1:(K+1)÷2 , j = 1:(J+1)÷2
                if j == (J+1)÷2
                    @test 2/sqrt(2) * M₁[2K-2k+1,2J-2j+1]M₁[2k-1,2j-1] ≈ M[Block(K,J)[k,j]]
                else
                    @test M₁[2K-2k+1,2J-2j+1]M₁[2k-1,2j-1] ≈ M[Block(K,J)[k,j]]
                end
            end
        end


        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (6,0)
        @test M₁[7,7]M₁[1,1] ≈ M[Block(4,4)[1,1]]
        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (4,2)
        @test M₁[7,5]M₁[1,3] ≈ M[Block(4,4)[1,2]]
        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (8,0)
        @test M₁[7,9]M₁[1,1] ≈ M[Block(4,5)[1,1]]
        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (6,2)
        @test M₁[7,7]M₁[1,3] ≈ M[Block(4,5)[1,2]]
        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (8,2)
        @test M₁[7,9]M₁[1,3] ≈ M[Block(4,6)[1,2]]


        # (2k,2j) = (4,2)
        # (2k̃,2j̃) = (4,2)
        @test M₁[5,5]M₁[3,3] + M₁[5,3]^2 ≈ M[Block(4,4)[2,2]]
        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (4,2)
        @test M₁[7,5]M₁[1,3] ≈ M[Block(4,4)[1,2]]
        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (8,0)
        @test M₁[7,9]M₁[1,1] ≈ M[Block(4,5)[1,1]]
        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (6,2)
        @test M₁[7,7]M₁[1,3] ≈ M[Block(4,5)[1,2]]
        # (2k,2j) = (6,0)
        # (2k̃,2j̃) = (8,2)
        @test M₁[7,9]M₁[1,3] ≈ M[Block(4,6)[1,2]]


        K,J = 6,6
        k,j = 1
        for j = 1:size(M[Block(K,J)],2)
            @test M₁[2K-1,2J-2j+1]M₁[2k-1,2j-1] ≈ M[Block(K,J)[k,j]]
        end


        Δ = R * weaklaplacian(W²)[Block.(1:2:2N), Block.(1:2:2N)] * R'
        D = weaklaplacian(W)

        @test 2M₁[1,1]D[1,1] ≈ Δ[1,1]
        @test  D[1,1]M₁[1,3]*sqrt(2) ≈ Δ[1,2]

        @test (M₁[3,3]D[1,1] + M₁[1,1]D[3,3]) ≈ Δ[2,2]
        @test M₁[3,5]D[1,1] ≈ Δ[2,3]
        @test 2/sqrt(2) * M₁[1,3]D[3,3] ≈ Δ[2,4]

        (K,J,k,j) = (3,3,1,1)
        # (2k,2j) = (4,0)
        # (2k̃,2j̃) = (4,0)
        @test M₁[5,5]D[1,1] + D[5,5]M₁[1,1] ≈ Δ[Block(K,J)[k,j]]
        
        
    
        for K = 1:2:9, J = 1:2:9
            for k=1:(K+1)÷2 , j = 1:(J+1)÷2
                if (j == (J+1)÷2 && k == (K+1)÷2)
                    @test 1/2 * (M₁[2K-2k+1,2J-2j+1]D[2k-1,2j-1] + M₁[2K-2k+1,2j-1]D[2J-2j+1,2k-1] + M₁[2k-1,2J-2j+1]D[2K-2k+1,2j-1] + M₁[2k-1,2j-1]D[2K-2k+1,2J-2j+1]) ≈ Δ[Block(K,J)[k,j]]
                else
                    @test (M₁[2K-2k+1,2J-2j+1]D[2k-1,2j-1] + M₁[2K-2k+1,2j-1]D[2J-2j+1,2k-1] + M₁[2k-1,2J-2j+1]D[2K-2k+1,2j-1] + M₁[2k-1,2j-1]D[2K-2k+1,2J-2j+1]) ≈ Δ[Block(K,J)[k,j]]
                end
            end
        end

        for K = 1:2:9, J = 2:2:9
            for k=1:(K+1)÷2 , j = 1:(J+1)÷2
                if (j == (J+1)÷2 && k == (K+1)÷2)
                    @test sqrt(2)/2 * (M₁[2K-2k+1,2J-2j+1]D[2k-1,2j-1] + M₁[2K-2k+1,2j-1]D[2J-2j+1,2k-1] + M₁[2k-1,2J-2j+1]D[2K-2k+1,2j-1] + M₁[2k-1,2j-1]D[2K-2k+1,2J-2j+1]) ≈ Δ[Block(K,J)[k,j]]
                else
                    @test (M₁[2K-2k+1,2J-2j+1]D[2k-1,2j-1] + M₁[2K-2k+1,2j-1]D[2J-2j+1,2k-1] + M₁[2k-1,2J-2j+1]D[2K-2k+1,2j-1] + M₁[2k-1,2j-1]D[2K-2k+1,2J-2j+1]) ≈ Δ[Block(K,J)[k,j]]
                end
            end
        end

        for K = 2:2:11, J = 2:2:11
            for k=1:(K+1)÷2, j = 1:(J+1)÷2
                @test (M₁[2K-2k+1,2J-2j+1]D[2k-1,2j-1] + M₁[2K-2k+1,2j-1]D[2J-2j+1,2k-1] + M₁[2k-1,2J-2j+1]D[2K-2k+1,2j-1] + M₁[2k-1,2j-1]D[2K-2k+1,2J-2j+1]) ≈ Δ[Block(K,J)[k,j]]
            end
        end

        for K = 2:2:9, J = 1:2:9
            for k=1:(K+1)÷2 , j = 1:(J+1)÷2
                if (j == (J+1)÷2 && k == (K+1)÷2)
                    @test sqrt(2)/2 * (M₁[2K-2k+1,2J-2j+1]D[2k-1,2j-1] + M₁[2K-2k+1,2j-1]D[2J-2j+1,2k-1] + M₁[2k-1,2J-2j+1]D[2K-2k+1,2j-1] + M₁[2k-1,2j-1]D[2K-2k+1,2J-2j+1]) ≈ Δ[Block(K,J)[k,j]]
                else
                    @test (M₁[2K-2k+1,2J-2j+1]D[2k-1,2j-1] + M₁[2K-2k+1,2j-1]D[2J-2j+1,2k-1] + M₁[2k-1,2J-2j+1]D[2K-2k+1,2j-1] + M₁[2k-1,2j-1]D[2K-2k+1,2J-2j+1]) ≈ Δ[Block(K,J)[k,j]]
                end
            end
        end


        
    end
end