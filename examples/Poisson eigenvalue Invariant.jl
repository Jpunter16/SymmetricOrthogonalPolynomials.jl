using SymmetricOrthogonalPolynomials, SparseArrays, Arpack, LinearAlgebra, CairoMakie, Combinatorics
import Arpack:eigs

splus1(n::Int) = -1 / ((2n + 3) * (2n + 5))
s0(n::Int) = 2 / ((2n + 1) * (2n + 5))
sminus1(n::Int) = -1 / ((2n + 1) * (2n + 3))

# |Stab(μ)| = ∏ mult! (zeros count as slots)
stab(μ) = prod(factorial(count(==(v), μ)) for v in unique(μ))

# admissible (up,stay,down), up+stay+down = m: shift 2/0/-2 needs degree >= 0/0/2
comps(v, m) = [
    (up, stay, m - up - stay) for up = 0:m for
    stay = 0:(m-up) if (m - up - stay == 0 || v >= 2)
]

# ---------- conversion operator: m_α^W = Σ_μ S_{αμ} m_μ^{3/2} ----------
function conversionCoeffs(α::Vector{Int})
    vals = sort(unique(α); rev=true)
    mult = Dict(v => count(==(v), α) for v in vals)
    buckets = Dict{Vector{Int},Float64}()
    for combo in Iterators.product((comps(v, mult[v]) for v in vals)...)
        w = 1 / 1
        landing = Int[]
        for (u, (up, stay, down)) in enumerate(combo)
            v = vals[u]
            w *=
                splus1(v)^up * s0(v)^stay * sminus1(v)^down /
                    (factorial(up) * factorial(stay) * factorial(down))
            append!(landing, fill(v + 2, up))
            append!(landing, fill(v, stay))
            append!(landing, fill(v - 2, down))
        end
        μ = sort(landing; rev=true)
        buckets[μ] = get(buckets, μ, 0 / 1) + w
    end
    Dict(μ => stab(μ) * w for (μ, w) in buckets if !iszero(w))   # note: +, no diff slot
end

function ConvMatInvar(N::Int, dim::Int)
    parts = reduce(vcat, Partition_n_parts(i, dim) for i = 0:(N-1))
    dict_ind = makeDictInd(parts)
    n_len = length(parts)
    diag_11 = [2/(2*n+3) for n = 0:(N-1)]
    I = Int[]
    J = Int[]
    V = Float64[]
    for i in eachindex(parts)
        lapcof=conversionCoeffs(parts[i]) #name of lapcof is not representative
        for j in keys(lapcof)
            if sum(j)<=N-1
                w=lapcof[j]*length(multiset_permutations(j))
                for part_i in eachindex(j)
                    w *= diag_11[j[part_i]+1]
                end
                push!(I, i)
                push!(J, dict_ind[j])
                push!(V, w)
            end
        end
    end
    sparse(I, J, V, n_len, n_len)
end

dims=[3,5,6,7]

N=ones(Int,length(dims)).*30

fig=Figure()

index=[fig[1,1],fig[1,2],fig[2,1],fig[2,2]]

for i in eachindex(dims)
    ax=Axis(index[i]; xlabel="Coefficients", ylabel="Coefficient magnitude",yscale=log10, title="Approximation convergence dim=$(dims[i])")
    dim=dims[i]
    S=LaplacianMat(N[i], dim)
    M= ConvMatInvar(N[i], dim)

    F = eigen(Symmetric(Matrix(S)), Symmetric(Matrix(M)))
    k = argmin(abs.(F.values))
    λ1 = F.values[k]        
    u1 = F.vectors[:, k]
    u1=u1/norm(u1)
    aux=1
    u_deg=Float64[]
    for j= 1:(N[i]-1)
        len=length(Partition_n_parts(j,dim))
        val=norm(u1[aux:aux+len-1])
        if val !=0.0
            append!(u_deg,val)
        end
        aux=aux+len
    end
    exactlam = -dim*(pi^2/4)
    println("Error dim=$dim: $(abs(λ1-exactlam))")
    #ylims!(extrema(abs.(u_deg))[1]+2.2e-16,extrema(abs.(u_deg))[2])
    scatter!(ax,1:length(u_deg),abs.(u_deg),label="dim=$(dims[i])")
end

save("poisson eigenvalue Invariant.pdf",fig)
fig