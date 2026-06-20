using Test, SymmetricOrthogonalPolynomials, StaticArrays

P=Legendre()
X=SVector(0.1,0.2,0.3)

Ps=S3Invariant(P)

Ps[X,Block(1)[1]]

n=0
@test  Ps[X,Block(1)[1]]== 1.0

n=1
@test Ps[X,Block(2)[1]] ≈  0.34641016

n=2
@test Ps[X,Block(3)[1]] ≈  -0.74478185
@test Ps[X,Block(3)[2]] ≈  0.063508530

n=3
@test Ps[X,Block(4)[1]] ≈  -0.46765372
@test Ps[X,Block(4)[2]] ≈  -0.21555510
@test Ps[X,Block(4)[3]] ≈  0.0060000000

n=4
@test Ps[X,Block(5)[1]] ≈  0.37116405
@test Ps[X,Block(5)[2]] ≈  -0.12267861
@test Ps[X,Block(5)[3]] ≈  0.31813443
@test Ps[X,Block(5)[4]] ≈-0.028636573