using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials, CairoMakie
import CairoMakie:spy!

N=20

alltups = vcat(([part for n in 0:N for part in Partition_2_parts(n)]))
inds_keep = findall(t -> t[1] >= 2 && t[2] >= 2, alltups)
alltups_keep = alltups[inds_keep]
inds_o = findall(t -> isodd(sum(t)), alltups)
inds_e = findall(t -> iseven(sum(t)), alltups)
inds=[inds_o; inds_e]

fig = Figure(size=(800,600))
Axis(fig[1,1]; yreversed=true, title="Not using symmetry")

Δ=(getLaplacianS2InvariantBasis(S2Invariant((Ultraspherical(-0.5))),N))
Δ[abs.(Δ) .< 1e-10] .= 0

spy!(Δ)

Axis(fig[1,2]; yreversed=true, title="Odd-even symmetry")


spy!(Δ[inds, inds])

Axis(fig[2,1]; yreversed=true, title="n=0,1 excluded ")

Δ_sym = Δ[inds_keep, inds_keep]
Δ_sym[abs.(Δ_sym) .< 1e-10] .= 0
spy!(Δ_sym)

Axis(fig[2,2]; yreversed=true, title="n=0,1 excluded, odd-even symmetry")

inds_combined = filter(∈(inds_keep), inds)
spy!(Δ[inds_combined, inds_combined])



str= "Sparsity of Laplace operator n=" * string(N) *".png"


save(str,fig)


fig