using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials, CairoMakie
import CairoMakie:spy!

N=20

inds_keep,inds_decouple=S3AntiAndInvariantBasisPruner(N)[1]
inds_keep_parity,inds_keep_Invariance_parity=S3AntiAndInvariantBasisPruner(N)[2]

fig = Figure(size=(800,800))
Axis(fig[1,1]; yreversed=true, title="Permutation and negation")

Δ=(getLaplacianS3AntiAndInvariantBasis(S3Invariant((Ultraspherical(-0.5))),N))

Δ[abs.(Δ) .< 1e-16] .= 0

spy!(Δ[inds_keep...,inds_keep...])
Axis(fig[1,2]; yreversed=true, title="Invariant and Anti-Invariant adapted")

spy!(Δ[inds_decouple...,inds_decouple...])

Axis(fig[2,1];yreversed=true,title="Negation adapted")
spy!(Δ[inds_keep_parity...,inds_keep_parity...])

Axis(fig[2,2]; yreversed=true, title="Invariant, Anti-Invariant and Negation adapted")

spy!(Δ[inds_keep_Invariance_parity...,inds_keep_Invariance_parity...])

str= "Sparsity of Laplace operator under S3 Invariant and Anti-invariant n=" * string(N) *".png"


save(str,fig)


fig