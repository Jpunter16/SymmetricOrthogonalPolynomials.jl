########
# Invariant polynomials with respect to S3 are given by
# P_{i}(x) P_{j}(y)P_{k}(z) +  "all other permuations of (x,y,z)"
# The first few are 
# 1
# -----
# 2P_1(x) + 2P_1(y)+ 2P_1(z)
# ----
# 2P_2(x) + 2P_2(y) + 2P_2(z)
# 2P_1(x)P_1(y) + 2P_1(x)P_1(z) + 2P_1(y)P_1(z)
# ----
# 2P_3(x) + 2P_3(y) + 2P_3(z)
# P_2(x)P_1(y)+P_2(x)P_1(z)+P_2(y)P_1(z) +P_2(z)P_1(y)+P_2(y)P_1(x) +P_2(z)P_1(x)
# 6P_1(x)P_1(y)P_1(z)
########

struct S3Invariant{T,B} <: MultivariateOrthogonalPolynomial{3,T}
    basis::B
end
S3Invariant(B::AbstractQuasiMatrix{T}) where T = S3Invariant{T, typeof(B)}(B)

S3Invariant() = S3Invariant(Legendre())


S3axis(K) = BlockedOneTo(round.((4:(K+4)) .^ 2 ./ 12))
axes(::S3Invariant) = (Inclusion(ChebyshevInterval() × ChebyshevInterval() × ChebyshevInterval()), S3axis(∞))




function getindex(Q::S3Invariant, 𝐱::SVector{3}, Kk::BlockIndex{1})
    x,y,z = 𝐱
    K,k = block(Kk), blockindex(Kk)
    K_int=Int(K)-1
    part= Partition_3_parts(K_int)[k]
    sigmaall = unique(collect(permutations(part)))
    numperms = length(sigmaall)
    norm = 1 / sqrt(numperms)
    
    total = 0
    P=Q.basis

    for i = 1:numperms
        mult=P[x,sigmaall[i][1]+1]*P[y,sigmaall[i][2]+1]*P[z,sigmaall[i][3]+1]
        total = total + mult
    end
    norm * total

    #=num_diff=(ℓ==μ)+(ℓ==ρ)+(μ==ρ)

    (P[x,ℓ+1]*P[y,μ+1]*P[z, ρ+1]+
    P[x,ℓ+1]*P[z,μ+1]*P[y, ρ+1]+
    P[z,ℓ+1]*P[y,μ+1]*P[x, ρ+1]+
    P[y,ℓ+1]*P[x,μ+1]*P[z, ρ+1]+
    P[z,ℓ+1]*P[x,μ+1]*P[y, ρ+1]+
    P[y,ℓ+1]*P[z,μ+1]*P[x, ρ+1])/sqrt(6+3*num_diff^2+num_diff) =#

    #Ask about normalization

    #(Q.basis[x,ℓ+1]Q.basis[y,μ+1]+Q.basis[x,μ+1]Q.basis[y,ℓ+1])/sqrt(2 + 2*(ℓ == μ)) # scaling is to ensure unitary change-of-basis
end

getindex(Q::S3Invariant, 𝐱::SVector{3}, k::Int) = Q[𝐱,findblockindex(axes(Q,2),k)]

getindex(Q::S3Invariant, 𝐱::SVector{3}, J::Block{1}) = [Q[𝐱,J[j]] for j = 1:length(axes(Q,2)[J])]
getindex(Q::S3Invariant, 𝐱::SVector{3}, JR::BlockOneTo) = mortar([Q[𝐱,J] for J in JR])


struct Laplacian{T} <: AbstractBandedBlockBandedMatrix{T}
    D::AbstractMatrix{T} # 1D Weak Laplacian
end

#=
struct S3KronVector{T,D<:AbstractVector{T}} <: AbstractBlockVector{T}
    d::D
end

axes(::S3KronVector) = (S3axis(∞),)
size(::S3KronVector) = (ℵ₀,)

function getindex(D::S3KronVector, K::Block{1})
    K̃ = Int(K)
    D.d[1:2:K̃] .* D.d[2K̃-1:-2:K̃]    
end
getindex(D::S3KronVector, Kk::BlockIndex{1}) = D[block(Kk)][blockindex(Kk)]
getindex(D::S3KronVector, k::Int) = D[findblockindex(axes(D,1), k)]




function grammatrix(Q::S3Invariant)
    M = grammatrix(Q.basis)
    Diagonal(S3KronVector(M.diag))
end

@simplify function *(Ac::QuasiAdjoint{<:Any,<:S3Invariant}, B::S3Invariant)
    M = (Ac').basis'B.basis
    Diagonal(S3KronVector(M.diag))
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

=#