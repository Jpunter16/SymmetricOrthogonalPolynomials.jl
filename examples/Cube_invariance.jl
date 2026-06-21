using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials, LazyBandedMatrices,CairoMakie, BlockDiagonals
import SparseArrays: sparse
import CairoMakie:spy!
import SymmetricOrthogonalPolynomials:cubeperm_inds


P3=S3Invariant(Ultraspherical(-0.5))

P=P3.basis


n=15  #degree of truncation +1 
N=n

P_diff=diff(P)

# 1D mass and stifness matrices
M1D = (P' * P)[1:N, 1:N]       #later change to avoid degree 0 and 1 plynomials that are not 0 at the boundary
S1D = (P_diff' * P_diff)[1:N, 1:N]    

# 3D stiffness matrix
∇ = sparse(KronTrav(S1D, M1D,M1D)) + sparse(KronTrav(M1D, S1D,M1D)+sparse(KronTrav(M1D,M1D,S1D)))

#jacobi matrix
J = jacobimatrix(P)[1:N, 1:N]
X = sparse(KronTrav(J,M1D,M1D))
Y = sparse(KronTrav(M1D,J,M1D))
Z = sparse(KronTrav(M1D,M1D,J))

V = (X-Y)^2 + (Y-Z)^2 + (X-Z)^2

# (x-y)^2 + (y-z)^2 + (x-z)^2

L = ∇ + V

Q=get_Q(N)
inds = cubeperm_inds(N)

alltups = vcat((lextuples(n) for n=1:N)...)
inds_o = findall(isodd,map(sum, alltups)); inds_e = findall(iseven,map(sum, alltups))
length(inds_o) # 95
length(inds_e) # 125
inds_eo = [inds_o; inds_e]

#ask about cubeperm_inds in original file (cubeperm.jl)
fig = Figure(size=(600,300))
Axis(fig[1,1]; yreversed=true, title="Reflection adapted")
spy!(sparse(Matrix((L)[inds_eo, inds_eo])))
Axis(fig[1,2]; yreversed=true, title="Permutation adapted")
spy!(round.(Matrix(Q'*L*Q)[inds,inds];digits=10))

str= "Sparsity of Schordinger operator n=" * string(n-1) *".png"

save(str,fig)

fig 
