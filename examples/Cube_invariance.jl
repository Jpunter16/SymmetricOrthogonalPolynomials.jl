using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials, LazyBandedMatrices,CairoMakie, BlockDiagonals
import SparseArrays: sparse
import CairoMakie:spy!

P3=S3Invariant(Ultraspherical(-0.5))

P=P3.basis


n=10  #degree of truncation +1 
#N=sum([binomial(i+2,3) for i=1:n])
N=n

P_diff=diff(P)

# 1D mass and stifness matrices
M1D = (P' * P)[1:N, 1:N]       
S1D = (P_diff' * P_diff)[1:N, 1:N]    

# 3D stiffness matrix
∇ = sparse(KronTrav(S1D, M1D,M1D)) + sparse(KronTrav(M1D, S1D,M1D)+sparse(KronTrav(M1D,M1D,S1D)))

#jacobi matrix
J=jacobimatrix(P)[1:N, 1:N]
X = sparse(KronTrav(J,M1D,M1D))
Y = sparse(KronTrav(M1D,J,M1D))
Z = sparse(KronTrav(M1D,M1D,J))

V = (X-Y)^2 + (Y-Z)^2 + (X-Z)^2

# (x-y)^2 + (y-z)^2 + (x-z)^2

L = ∇ + V

Q=get_Q(N)

#ask about cubeperm_inds in original file (cubeperm.jl)
fig = Figure(size=(600,300))
Axis(fig[1,1]; yreversed=true, title="Reflection adapted")
spy!(sparse(Matrix((L))))
Axis(fig[1,2]; yreversed=true, title="Permutation adapted")
spy!(round.(Matrix(Q'*L*Q);digits=10))
fig 


