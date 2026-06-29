using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials, CairoMakie
import CairoMakie:spy!

N=5

alltups = vcat((Partition_2_parts(n) for n=0:N))
inds_o = findall(isodd,map(sum, alltups)); inds_e = findall(iseven,map(sum, alltups))
inds_eo = [inds_o; inds_e]

fig = Figure(size=(600,300))
Axis(fig[1,1]; yreversed=true)

Δ=Matrix(getLaplacianS2InvariantBasis(S2Invariant((Ultraspherical(0.5))),N))

spy!(Δ)

Axis(fig[1,2]; yreversed=true)


spy!(Δ[inds_eo, inds_eo])
fig