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


S3axis(K) = blockedrange(((3:(K+2)) .^ 2 .+ 6) .÷ 12)
axes(::S3Invariant) = (Inclusion(ChebyshevInterval() × ChebyshevInterval() × ChebyshevInterval()), S3axis(∞))

function getexpresionS3Invariant(Q::S3Invariant,part::Partition3) #only gets polinomials of type P[x,n]*P[x,m], only taking the same space variable, not useful
    P=Q.basis
    sigmaall = unique(collect(permutations(part.p)))
    numperms = length(sigmaall)
    norm = 1 / sqrt(numperms)

    total = P[:,sigmaall[1][1]+1].*P[:,sigmaall[1][2]+1].*P[:,sigmaall[1][3]+1]
    for i = 2:numperms
        mult=P[:,sigmaall[i][1]+1].*P[:,sigmaall[i][2]+1].*P[:,sigmaall[i][3]+1]
        total = total .+ mult
    end
    norm * total
end

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

function checkPosiblePartitionLaplacian(v::Partition3)
    D= [zeros(Int,3) for _ in 1:3]
    for i=1:3
        aux=v.p
        aux[i]=aux[i]-2
        aux=sort(aux, rev=true)
        D[i]=aux
    end
    [all(x ->x >=0, v) for v in D]
end



function getLaplacianFromUltraspherical(Q::S3Invariant, n::Int) #get krontrav matrices
    D2=S3axis(n)
    C=Q.basis
    T = typeof(Q).parameters[1]
    m=2
    #this code from ultraspherical.jl line 152, differentiate ultraspherical
    μ = pochhammer(convert(T,C.λ),m)*convert(T,2)^m
    D2_mult = _BandedMatrix(Fill(μ,1,∞), ℵ₀, -m, m)
    Px=diff(C,2)
    Py=Px
    Pz=Px
    for i=0:n
        part=Partition_3_parts(i)
        for j=1:eachindex(part)
            part_lap=checkPosiblePartitionLaplacian(part[j])
            for k=1:eachindex(part_lap)
                if part_lap[k]
                end
            end
        end
    end
    D2_mult

end



struct S3InvariantLaplacian{T,B} <: MultivariateOrthogonalPolynomial{3,T}
    basis::B
    D::BlockArray{T} # Laplacian
end

function getindex(Δ::S3InvariantLaplacian, X::SVector{3}, Kk::BlockIndex{1})

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