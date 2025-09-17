using LinearAlgebra, ClassicalOrthogonalPolynomials, Test


###
# Legendre
###

P = Legendre()
𝕍 = (n,x,y) -> hcat([[P[x,n-k+1]*P[y,k+1], 0] for k=0:n]..., [[0,P[x,n-k+1]*P[y,k+1]] for k=0:n]...)

x,y=0.1,0.2;
n = 0
@test [0 1; 1 0] * 𝕍(0,y,x) ≈ 𝕍(0,x,y) * [0 1; 1 0]
@test [1 0; 0 -1] * 𝕍(0,x,-y) ≈ 𝕍(0,x,y) * [1 0; 0 -1]

n = 1
@test [0 1; 1 0] * 𝕍(n,y,x) ≈ 𝕍(n,x,y) * I(2(n+1))[end:-1:1,:]
@test [1 0; 0 -1] * 𝕍(n,x,-y) ≈ 𝕍(n,x,y) * Diagonal([(-1).^(0:n); -(-1).^(0:n)])



σ = I(2(n+1))[end:-1:1,:]
μ = Diagonal([(-1).^(0:n); -(-1).^(0:n)])

# σ * Q 

nullspace(Matrix([σ - I;  μ - I]))
nullspace(Matrix([σ - I;  μ + I]))
nullspace(Matrix([σ + I;  μ - I]))
nullspace(Matrix([σ + I;  μ + I]))

ℚ = (n,x,y) -> 𝕍(n,x,y) * [1 0 -1  0;
                           0 1  0  1; 
                           0 1  0 -1; 
                           1 0  1  0]/sqrt(2)

@test [0 1; 1 0] * ℚ(n,y,x) ≈ ℚ(n,x,y) * Diagonal([1,1,-1,-1])
@test [1 0; 0 -1] * ℚ(n,x,-y) ≈ ℚ(n,x,y) * Diagonal([1,-1,1,-1])


n = 2

σ = I(2(n+1))[end:-1:1,:]
μ = Diagonal([(-1).^(0:n); -(-1).^(0:n)])

@test [0 1; 1 0] * 𝕍(n,y,x) ≈ 𝕍(n,x,y) * σ
@test [1 0; 0 -1] * 𝕍(n,x,-y) ≈ 𝕍(n,x,y) * μ

m = size(σ,2)

Q = nullspace(Matrix([kron(I(2),σ) - kron([1 0; 0 -1],I(m));  kron(I(2),μ) - kron([0 1; 1 0],I(m))]))
Q = reshape(vec(Q), m, m)*sqrt(2)
@test Q'σ*Q ≈ blockdiag(fill(sparse([1 0; 0 -1]), m ÷2)...)
@test Q'μ*Q ≈ blockdiag(fill(sparse([0 1; 1 0]), m ÷2)...)

###
# Different combo
###

P = Legendre()
𝕍 = (n,x,y) -> hcat([P[x,n-k+1]*P[y,k+1]*[1, 1] for k=0:n]..., [P[x,n-k+1]*P[y,k+1]*[1,-1] for k=0:n]...)

x,y=0.1,0.2;
n = 0
@test [1 0; 0 -1] * 𝕍(0,x,-y) ≈ 𝕍(0,x,y) * [0 1; 1 0]
@test [0 1; 1 0] * 𝕍(0,y,x) ≈ 𝕍(0,x,y) * [1 0; 0 -1]

n = 1

@test [1 0; 0 -1] * 𝕍(n,x,-y) ≈ 𝕍(n,x,y) * [0 0 1 0
                                            0 0 0 -1
                                            1 0 0 0
                                            0 -1 0 0]
@test [0 1; 1 0] * 𝕍(n,y,x) ≈ 𝕍(n,x,y) *   [0 1 0 0
                                            1 0 0 0
                                            0 0 0 -1
                                            0 0 -1 0]






σ = [0 0 1 0
    0 0 0 -1
    1 0 0 0
    0 -1 0 0]
μ = [0 1 0 0
    1 0 0 0
    0 0 0 -1
    0 0 -1 0]


# σ * Q 

nullspace(Matrix([σ - I;  μ - I]))
nullspace(Matrix([σ - I;  μ + I]))
nullspace(Matrix([σ + I;  μ - I]))
nullspace(Matrix([σ + I;  μ + I]))

ℚ = (n,x,y) -> 𝕍(n,x,y) * [1 -1 -1  1;
                           1 -1  1  -1; 
                           1 1  -1 -1; 
                           -1 -1  -1  -1]/2

@test [0 1; 1 0] * ℚ(n,y,x) ≈ ℚ(n,x,y) * Diagonal([1,1,-1,-1])
@test [1 0; 0 -1] * ℚ(n,x,-y) ≈ ℚ(n,x,y) * Diagonal([1,-1,1,-1])


n = 2

σ = [0 0 1 0 0 0
    0 1 0 0 0 0
    1 0 0 0 0 0
    0 0 0 0 0 -1
    0 0 0 0 -1 0
    0 0 0 -1 0 0]
μ = [0 0  0 1 0 0
0 0  0 0 -1 0
0 0  0 0 0 1
1 0  0 0 0 0
0 -1 0 0 0 0
0 0  1 0 0 0]

@test [0 1; 1 0] * 𝕍(n,y,x) ≈ 𝕍(n,x,y) *   σ
@test [1 0; 0 -1] * 𝕍(n,x,-y) ≈ 𝕍(n,x,y) * μ

m = size(σ,2)

Q = nullspace(Matrix([kron(I(2),σ) - kron([1 0; 0 -1],I(m));  kron(I(2),μ) - kron([0 1; 1 0],I(m))]))
Q = reshape(vec(Q), m, m)*sqrt(2)
@test Q'σ*Q ≈ blockdiag(fill(sparse([1 0; 0 -1]), m ÷2)...)
@test Q'μ*Q ≈ blockdiag(fill(sparse([0 1; 1 0]), m ÷2)...)


###
# Mixed
###
W = Weighted(Jacobi(1,1))
𝕍 = (n,x,y) -> hcat([[P[x,n-k+1]*W[y,k+1], 0] for k=0:n]..., [[0,W[x,n-k+1]*P[y,k+1]] for k=0:n]...)

x,y=0.1,0.2;
n = 0
@test [0 1; 1 0] * 𝕍(0,y,x) ≈ 𝕍(0,x,y) * [0 1; 1 0]
@test [1 0; 0 -1] * 𝕍(0,x,-y) ≈ 𝕍(0,x,y) * [1 0; 0 -1]

n = 1
@test [0 1; 1 0] * 𝕍(n,y,x) ≈ 𝕍(n,x,y) * I(2(n+1))[end:-1:1,:]
@test [1 0; 0 -1] * 𝕍(n,x,-y) ≈ 𝕍(n,x,y) * Diagonal([(-1).^(0:n); -(-1).^(0:n)])



σ = I(2(n+1))[end:-1:1,:]
μ = Diagonal([(-1).^(0:n); -(-1).^(0:n)])

# σ * Q 

nullspace(Matrix([σ - I;  μ - I]))
nullspace(Matrix([σ - I;  μ + I]))
nullspace(Matrix([σ + I;  μ - I]))
nullspace(Matrix([σ + I;  μ + I]))

ℚ = (n,x,y) -> 𝕍(n,x,y) * [1 0 -1 0; 0 1 0 1; 0 1 0 -1; 1 0 1 0]/sqrt(2)

@test [0 1; 1 0] * ℚ(n,y,x) ≈ ℚ(n,x,y) * Diagonal([1,1,-1,-1])
@test [1 0; 0 -1] * ℚ(n,x,-y) ≈ ℚ(n,x,y) * Diagonal([1,-1,1,-1])


n = 2
@test [0 1; 1 0] * 𝕍(n,y,x) ≈ 𝕍(n,x,y) * I(2(n+1))[end:-1:1,:]
@test  [1 0; 0 -1] * 𝕍(n,x,-y) ≈ 𝕍(n,x,y) * Diagonal([(-1).^(0:n); -(-1).^(0:n)])

