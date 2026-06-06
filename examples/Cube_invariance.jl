using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials

P3=S3Invariant(Ultraspherical(-0.5))

P=P3.basis


N=10

Q=get_Q(N)

P_diff=diff(P)

# 1D mass and stifness matrices
M1D = (P' * P)[1:N, 1:N]       
S1D = (P_diff' * P_diff)[1:N, 1:N]    

# 3D stiffness matrix
S3D = sparse(KronTrav(S1D, M1D,M1D)) + sparse(KronTrav(M1D, S1D,M1D)+sparse(KronTrav(M1D,M1D,S1D)))



