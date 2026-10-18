using SymmetricOrthogonalPolynomials, SparseArrays, Arpack, LinearAlgebra, CairoMakie
import Arpack:eigs

function LaplacianMatAnti(N::Int, dim::Int)
    parts=reduce(vcat, Partition_diff_n_parts(i, dim) for i = 0:(N-1))
    dict_ind = makeDictInd(parts)
    n_len=length(parts)
    diag_11=[2/(2*n+3) for n = 0:(N-1)]
    I = Int[]
    J = Int[]
    V = Float64[]
    for i in eachindex(parts)
        lapcof=laplacianCoeffsAlter(parts[i])
        for j in keys(lapcof)
            if sum(j)<=N-1
                #w=factorial(dim)* lapcof[j]
                w=lapcof[j]
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

function ConvMatAnti(N::Int, dim::Int)
    parts=reduce(vcat, Partition_diff_n_parts(i, dim) for i = 0:(N-1))
    dict_ind = makeDictInd(parts)
    n_len=length(parts)
    diag_11=[2/(2*n+3) for n = 0:(N-1)]
    I = Int[]
    J = Int[]
    V = Float64[]
    for i in eachindex(parts)
        lapcof=convCoeffsAlter(parts[i])
        for j in keys(lapcof)
            if sum(j)<=N-1
                #w=factorial(dim)* lapcof[j]
                w=lapcof[j]
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

# Lowest eigenpair of the symmetric generalized problem S v = λ M v.
# RQI on its own converges to whichever eigenvalue is nearest the start, so first
# run shift-0 inverse iteration (→ eigenvalue of smallest |λ|) and then polish with RQI.
function lowest_eig(S, M; tol=1e-12, maxit=1000)
    v = zeros(size(S,1)); v[1] = 1
    v /= sqrt(v'*M*v)
    lam = v'*S*v
    F = lu(S)
    for k = 1:maxit
        w = F \ (M*v)
        v = w / sqrt(w'*M*w)
        lamnew = v'*S*v
        done = abs(lamnew-lam) < 1e-4*abs(lamnew)
        lam = lamnew
        done && break
    end
    for k = 1:maxit
        w = (S - lam*M) \ (M*v)
        v = w / sqrt(w'*M*w)
        lamnew = v'*S*v
        done = abs(lamnew-lam) < tol*abs(lamnew)
        lam = lamnew
        done && break
    end
    lam, v
end

f(n)=round(Int,n*(n-1)/2)

dims=[4,5,6,7]

N_true=31

N=map(f, dims).+N_true

fig=Figure()

mark=[:circle, :utriangle, :rect, :hexagon]

ax=Axis(fig[1,1]; xlabel="Degree", ylabel="Coefficient magnitude",yscale=log10, title="Magnitude of coefficients per total degree: Fermion case")

for i in eachindex(dims)
    dim=dims[i]
    S=sparse(LaplacianMatAnti(N[i], dim)')
    M=sparse(ConvMatAnti(N[i], dim)')

    λ1, u1 = lowest_eig(S, M)
    u1/=norm(u1)
    aux=1
    degs=round(Int,dim*(dim-1)/2):(N[i]-1)
    u_deg=Float64[]
    for j in degs
        len=length(Partition_diff_n_parts(j,dim))
        append!(u_deg,norm(u1[aux:aux+len-1]))
        aux=aux+len
    end
    exactlam = -(pi^2/4)*sum([i^2 for i in 1:dim])
    println("Error dim=$dim: $(abs(λ1-exactlam))")
    #ylims!(extrema(abs.(u_deg))[1]+2.2e-16,extrema(abs.(u_deg))[2])
    # keep every degree, but only plot the non-negligible ones, at their true position
    keep=u_deg .> 1e-15
    xs=(0:length(u_deg)-1)[keep]
    scatterlines!(ax,xs,u_deg[keep]; marker=mark[i], label="dim=$(dims[i])", strokewidth=1, strokecolor=:black)
end
axislegend(ax, position=:lb)
save("poisson eigenvalue anti-invariant.pdf",fig)
fig