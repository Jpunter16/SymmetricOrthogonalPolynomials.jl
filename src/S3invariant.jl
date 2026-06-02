########
# Invariant polynomials with respect to S3 are given by
# P_{i}(x) P_{j}(y)P_{k}(z) +  "all other permuations of (x,y,z)"
# The first few are 
# 1
# -----
# P_1(x) + P_1(y)+P_1(z)
# ----
# P_2(x) + P_2(y)+P_2(z)
# 2P_1(x)P_1(y)+2P_1(x)P_1(z)+2P_1(y)P_1(z)
########

struct S3Invariant{T,B} <: MultivariateOrthogonalPolynomial{3,T}
    basis::B
end
S3Invariant(B::AbstractQuasiMatrix{T}) where T = S3Invariant{T, typeof(B)}(B)

S3Invariant() = S3Invariant(Normalized(Legendre()))


#continue here
S3axis(K) = BlockedOneTo((2:(K+1)) .^ 2 .÷ 4)
axes(::S3Invariant) = (Inclusion(ChebyshevInterval() × ChebyshevInterval()), dihedralaxis(∞))




function getindex(Q::S3Invariant, 𝐱::SVector{2}, Kk::BlockIndex{1})
    x,y = 𝐱
    K,k = block(Kk), blockindex(Kk)
    ℓ = 2*(Int(K)-k)
    μ = 2*(k-1)
    (Q.basis[x,ℓ+1]Q.basis[y,μ+1]+Q.basis[x,μ+1]Q.basis[y,ℓ+1])/sqrt(2 + 2*(ℓ == μ)) # scaling is to ensure unitary change-of-basis
end

getindex(Q::S3Invariant, 𝐱::SVector{2}, k::Int) = Q[𝐱,findblockindex(axes(Q,2),k)]

getindex(Q::S3Invariant, 𝐱::SVector{2}, J::Block{1}) = [Q[𝐱,J[j]] for j = 1:length(axes(Q,2)[J])]
getindex(Q::S3Invariant, 𝐱::SVector{2}, JR::BlockOneTo) = mortar([Q[𝐱,J] for J in JR])

struct DihedralKronVector{T,D<:AbstractVector{T}} <: AbstractBlockVector{T}
    d::D
end

axes(::DihedralKronVector) = (dihedralaxis(∞),)
size(::DihedralKronVector) = (ℵ₀,)


function getindex(D::DihedralKronVector, K::Block{1})
    K̃ = Int(K)
    D.d[1:2:K̃] .* D.d[2K̃-1:-2:K̃]    
end
getindex(D::DihedralKronVector, Kk::BlockIndex{1}) = D[block(Kk)][blockindex(Kk)]
getindex(D::DihedralKronVector, k::Int) = D[findblockindex(axes(D,1), k)]




function grammatrix(Q::S3Invariant)
    M = grammatrix(Q.basis)
    Diagonal(DihedralKronVector(M.diag))
end

@simplify function *(Ac::QuasiAdjoint{<:Any,<:S3Invariant}, B::S3Invariant)
    M = (Ac').basis'B.basis
    Diagonal(DihedralKronVector(M.diag))
end



function dihedralconversion(N)
    R = BlockBandedMatrix{Float64}(undef, (dihedralaxis(N), BlockedOneTo(cumsum(1:2:2N))), (0,0)); fill!(R, 0)
    for K = 1:N
        for k = 1:(K÷2)
            R[Block(K,K)[k,2k-1]] = R[Block(K,K)[k,2K-2k+1]] = 1/sqrt(2)
        end
        isodd(K) && (R[Block(K,K)[K÷2+1,K]] = 1)
    end
    R
end

struct DihedralWeakLaplacian{T} <: AbstractBandedBlockBandedMatrix{T}
    D::AbstractMatrix{T} # 1D Weak Laplacian
    M::AbstractMatrix{T} # 1D Mass Matrix
end


axes(::DihedralWeakLaplacian) = (dihedralaxis(∞),dihedralaxis(∞))

blockbandwidths(::DihedralWeakLaplacian) = (1,1)
subblockbandwidths(::DihedralWeakLaplacian) = (1,1)

function viewblock(Δ::DihedralWeakLaplacian, KJ::Block{2})
    D,M = Δ.D,Δ.M
    K,J = KJ.n
    m,n = (K+1)÷2, (K+1)÷2
    if K == J
        BandedMatrix(0 => [M[2K-2k+1,2K-2k+1]D[2k-1,2k-1] + 2M[2K-2k+1,2k-1]D[2K-2k+1,2k-1] + M[2k-1,2k-1]D[2K-2k+1,2K-2k+1] for k=1:m])
    elseif J == K+1

    elseif K == J+1

    else

    end
end