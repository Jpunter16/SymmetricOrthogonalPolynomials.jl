using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials, CairoMakie
import CairoMakie:spy!

N=20

alltups = vcat(([part.p for n in 0:N for part in Partition_2_parts(n)]))
inds_o = findall(isodd,map(sum, alltups)); inds_e = findall(iseven,map(sum, alltups))
inds_eo = [inds_o; inds_e]

fig = Figure(size=(800,400))
Axis(fig[1,1]; yreversed=true, title="Not using symmetry")

Δ=Matrix(getLaplacianS2InvariantBasis(S2Invariant((Ultraspherical(0.5))),N))

spy!(Δ)

Axis(fig[1,2]; yreversed=true, title="Odd even symmetry")


spy!(Δ[inds_eo, inds_eo])


str= "Sparsity of Laplace operator n=" * string(N) *".png"


save(str,fig)


fig