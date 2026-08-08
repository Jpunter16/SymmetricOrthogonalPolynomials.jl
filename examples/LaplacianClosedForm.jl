using SymmetricOrthogonalPolynomials, ClassicalOrthogonalPolynomials, CairoMakie

delt3=getLaplacianClosedForm((Ultraspherical(-0.5)[:,3:end]),20,3)

fig = Figure(size=(800,800))
Axis(fig[1,1]; yreversed=true, title="Permutation of S3")

spy!(delt3)


delt4=getLaplacianClosedForm((Ultraspherical(-0.5)[:,3:end]),20,4)

Axis(fig[1,2];yreversed=true, title="Permutation of S4")


spy!(delt4)


delt5=getLaplacianClosedForm((Ultraspherical(-0.5)[:,3:end]),20,5)

Axis(fig[2,1];yreversed=true, title="Permutation of S5")

spy!(delt5)

delt6=getLaplacianClosedForm((Ultraspherical(-0.5)[:,3:end]),20,6)

Axis(fig[2,2];yreversed=true, title="Permutation of S6")

spy!(delt6)

save("Closed Form S3,S4,S5,S6 n=20.png",fig)

fig