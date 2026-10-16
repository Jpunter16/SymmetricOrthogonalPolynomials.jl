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

f(n)=round(Int,n*(n-1)/2)

dims=[3,5,6,7]

N=map(f, dims).+25

fig=Figure()

index=[fig[1,1],fig[1,2],fig[2,1],fig[2,2]]

for i in eachindex(dims)
    ax=Axis(index[i]; xlabel="Coefficients", ylabel="Coefficient magnitude",yscale=log10, title="Approximation convergence dim=$(dims[i])")
    dim=dims[i]
    S=LaplacianMatAnti(N[i], dim)
    M= ConvMatAnti(N[i], dim)

    F = eigen(Symmetric(Matrix(S)), Symmetric(Matrix(M)))
    k = argmin(abs.(F.values))
    λ1 = F.values[k]        
    u1 = F.vectors[:, k]
    u1=u1/norm(u1)
    aux=1
    u_deg=Float64[]
    for j= round(Int,dim*(dim-1)/2):(N[i]-1)
        len=length(Partition_diff_n_parts(j,dim))
        val=norm(u1[aux:aux+len-1])
        if val !=0.0
            append!(u_deg,val)
        end
        aux=aux+len
    end
    exactlam = -(pi^2/4)*sum([i^2 for i in 1:dim])
    println("Error dim=$dim: $(abs(λ1-exactlam))")
    #ylims!(extrema(abs.(u_deg))[1]+2.2e-16,extrema(abs.(u_deg))[2])
    scatter!(ax,1:length(u_deg),abs.(u_deg),label="dim=$(dims[i])")
end

save("poisson eigenvalue anti-invariant.pdf",fig)
fig