using LinearAlgebra, DynamicPolynomials, NumericalRepresentationTheory, LazyBandedMatrices
import SymmetricOrthogonalPolynomials: cubevectorrepresentation, cubevector_filter


@polyvar(x,y,z)


Q = Array{Matrix{Float64}}(undef, 10)


ℙ₂ = (n,x,y) -> [legendrep(n-k,x)legendrep(k,y) for k=0:n]
ℙ = (n,x,y,z) -> vcat([legendrep(n-k,x) .* ℙ₂(k,y,z) for k=0:n]...)

𝕍 = (n,x,y,z) -> hcat(vcat.(ℙ(n,x,y,z), 0,0)..., vcat.(0, ℙ(n,x,y,z),0)...,vcat.(0, 0, ℙ(n,x,y,z))...)

𝕍(2,x,y,z)


s₁ = [0 0 1;
0 -1 0;
1 0 0]

s₂ = [-1 0 0;
0 0 -1;
0 -1 0]

s₃ = [0 0 -1;
0 -1 0;
-1 0 0]

for n = 1:10
    ρ = cubevectorrepresentation(n)
    S₁,S₂,S₃ = ρ.generators


    @test s₁*𝕍(n-1,s₁*[x,y,z]...) == 𝕍(n-1,x,y,z) * Matrix(S₁')
    @test s₂*𝕍(n-1,s₂*[x,y,z]...) == 𝕍(n-1,x,y,z) * Matrix(S₂')
    @test s₃*𝕍(n-1,s₃*[x,y,z]...) == 𝕍(n-1,x,y,z) * Matrix(S₃')


    @test S₁^2 == S₂^2 == S₃^2 == I
    @test S₁*S₃ == S₃*S₁
    @test (S₁*S₂)^3 == (S₂*S₃)^3 == I
end


N = 10
Qs = Array{Matrix{Float64}}(undef, N)
for n = 1:N
    @show n
    _,Qs[n] = blockdiagonalize(cubevectorrepresentation(n))
end


Q = blockdiag(sparse.(Qs)...)
P = Legendre()
D¹ = (P'diff(P))[1:N,1:N]
D² = (diff(P)'diff(P))[1:N,1:N]
M = (P'P)[1:N,1:N]


D_xy = sparse(KronTrav(M, D¹, D¹'))
D_xz = sparse(KronTrav(D¹, M, D¹'))
D_yz = sparse(KronTrav(D¹, D¹', M))
D_xx = sparse(KronTrav(M, M, D²))
D_yy = sparse(KronTrav(M, D², M))
D_zz = sparse(KronTrav(D², M, M))


Ng = 20
𝐱,𝐲,𝐳 = ClassicalOrthogonalPolynomials.grid(P, (Ng,Ng,Ng))
Pl = plan_transform(P, (Ng,Ng,Ng))
f = (x,y,z) -> 1 + x + y + z + x^2 + x*y + y^2 + x^3 + x^2*y + x*y^2 + y^3 + x^4 + z^4
g = (x,y,z) -> 1 + 2x + 3y + 4z + 4x^2 + 5x*y + 6y^2 + 7x^3 + x^2*y + x*y^2 + y^3 + y^4
h = (x,y,z) -> 1 + 2x - 3y + 4z + 4x^2 + 5x*y + 7y^2 - 7x^3 + x^2*y + x*y^2 + y^3 + y^4 + x*y*z^2
𝐟 = DiagTrav(Pl*f.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]
𝐠 = DiagTrav(Pl*g.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]
𝐡 = DiagTrav(Pl*h.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]

using ForwardDiff
f_y = (x,y,z) -> ForwardDiff.gradient(splat(f), [x,y,z])[2]
f_z = (x,y,z) -> ForwardDiff.gradient(splat(f), [x,y,z])[3]
g_x = (x,y,z) -> ForwardDiff.gradient(splat(g), [x,y,z])[1]
g_z = (x,y,z) -> ForwardDiff.gradient(splat(g), [x,y,z])[3]
h_x = (x,y,z) -> ForwardDiff.gradient(splat(h), [x,y,z])[1]
h_y = (x,y,z) -> ForwardDiff.gradient(splat(h), [x,y,z])[2]
h_z = (x,y,z) -> ForwardDiff.gradient(splat(h), [x,y,z])[3]

𝐟_y = DiagTrav(Pl*f_y.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]
𝐟_z = DiagTrav(Pl*f_z.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]
𝐠_x = DiagTrav(Pl*g_x.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]
𝐠_z = DiagTrav(Pl*g_z.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]
𝐡_x = DiagTrav(Pl*h_x.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]
𝐡_y = DiagTrav(Pl*h_y.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]
𝐡_z = DiagTrav(Pl*h_z.(𝐱,𝐲', reshape(𝐳,1,1,:)))[Block.(1:N)]

cubesum(f) = 8*(Pl * f.(𝐱,𝐲', reshape(𝐳,1,1,:)))[1]


@test 𝐟'D_yy*𝐟 ≈ cubesum((x,y,z) -> f_y(x,y,z)^2)
@test 𝐟'D_yy*𝐡 ≈ cubesum((x,y,z) -> f_y(x,y,z)h_y(x,y,z))
@test 𝐠'D_xx*𝐡 ≈ cubesum((x,y,z) -> g_x(x,y,z)h_x(x,y,z))
@test 𝐠'D_zz*𝐟 ≈ cubesum((x,y,z) -> g_z(x,y,z)f_z(x,y,z))
@test 𝐡'D_xy*𝐟 ≈ cubesum((x,y,z) -> h_x(x,y,z)f_y(x,y,z))
@test 𝐟'D_yz*𝐠 ≈ cubesum((x,y,z) -> f_y(x,y,z)g_z(x,y,z))
@test 𝐡'D_xz*𝐠 ≈ cubesum((x,y,z) -> h_x(x,y,z)g_z(x,y,z))


@test 𝐟'D_xx*𝐟 ≈ cubesum((x,y,z) -> f_x(x,y,z)^2)
cubesum((x,y,z) -> f(x,y,z)^2)

𝐟'KronTrav(M,M,M)[Block.(1:N),Block.(1:N)]*𝐟

L = [D_yy+D_zz  -D_xy'      -D_xz'; 
     -D_xy      D_xx+D_zz -D_yz'; 
     -D_xz     -D_yz      D_xx+D_yy]

@test L ≈ L'

@test cubesum((x,y,z) -> (h_y(x,y,z) - g_z(x,y,z))^2) + cubesum((x,y,z) -> (f_z(x,y,z) - h_x(x,y,z))^2) + cubesum((x,y,z) -> (g_x(x,y,z) - f_y(x,y,z))^2) ≈ [𝐟; 𝐠; 𝐡]'L*[𝐟; 𝐠; 𝐡]



ax = axes(KronTrav(M, M, D²),1)

# permute blocks to group by order
lngs = blocklengths(ax)
μ = length(lngs)
Σ = BlockArray{Int}(undef, 3*lngs, repeat(lngs,3)); Σ .= 0

for K = 1:μ
    m = lngs[K]
    Σ[Block(K,K)[1:m,1:m]] = Eye{Int}(m)
    Σ[Block(K,μ+K)[m+1:2m,1:m]] = Eye{Int}(m)
    Σ[Block(K,2μ+K)[2m+1:3m,1:m]] = Eye{Int}(m)
end

@test Σ'Σ == I




inds = Int[]
for s in (true,false), p in partitions(4)
    for j = 1:hooklength(p)
        append!(inds, findall(cubevector_filter((p,s), N, j)))
    end
end

sparse(round.(Q'*Σ*L*Σ'Q;digits=3)[inds,inds])

M³ = sparse(KronTrav(M,M,M))
sparse(round.(Q'*Σ*blockdiag(M³,M³,M³)*Σ'Q;digits=3))


BlockArray(round.(Q'*Σ*L*Σ'Q;digits=3), size.(Qs,1), size.(Qs,2))[Block(2,2)]


inds_t1 = findall(cubevector_filter((Partition(4), true), N))


L_p = (Q'*Σ*L*Σ'Q)[inds,inds]
M_p = (Q'*Σ*blockdiag(M³,M³,M³)*Σ'Q)[inds,inds]

sparse(round.(L_p+M_p;digits=3))[4:22,4:22] |> Matrix


eigen(L_p[4:22,4:22],M_p[4:22,4:22])

eigen(L_p[inds_t1,inds_t1], M_p[inds_t1,inds_t1])

eigvals(Symmetric(Matrix(L)), Symmetric(Matrix(blockdiag(M³,M³,M³))))


findall(≠(0),round.((L_p + M_p)[4,:];digits=3))

