using ClassicalOrthogonalPolynomials, Combinatorics, SymmetricOrthogonalPolynomials, SparseArrays, Plots, CairoMakie, Printf

function RHS(dim::Int)
    return X::AbstractVector -> begin
        @assert length(X)==dim "Expected lenght $dim"
        -dim*pi^2*prod([sinpi(x) for x in X])
    end 
end

function uExact(dim::Int)
    return X::AbstractVector -> begin
        @assert length(X)==dim "Expected lenght $dim"
        prod([sinpi(x) for x in X])
    end 
end

function assembleF(N,fhat::Dict{Vector{Int}, Float64}, dim)
    parts=reduce(vcat, Partition_n_parts(i, dim) for i=0:N-1)
    dict_ind=makeDictInd(parts)
    F=zeros((length(parts)))
    for i in keys(fhat)
        F[dict_ind[i]]=fhat[i]
    end
    F
end

function assembleF(N::Int,dim::Int)
    scale = -dim*pi^2
    Q = Ultraspherical(3/2)
    xg = ClassicalOrthogonalPolynomials.grid(Q[:, 1:N])                
    Vm = Q[xg, 1:N]                     # N×N generalized Vandermonde
    c = Vm \ sinpi.(xg)

    fhat = Dict{Vector{Int}, Float64}()

    for i in 0:N-1
        for μ in Partition_n_parts(i, dim)      
            fhat[μ] = scale * prod(c[k+1] for k in μ)
        end
    end
    assembleF(N,fhat,dim)
end

function makeDictInd(parts::Vector{Vector{Int}})
    d=Dict{Vector{Int},Int64}()
    for i in eachindex(parts)
        d[parts[i]]=i
    end
    d
end

function assembleA(N::Int, dim::Int)
    parts=reduce(vcat, Partition_n_parts(i, dim) for i=0:N-1)
    dict_ind=makeDictInd(parts)
    I=Int[]
    J=Int[]
    V=Float64[]
    for i in eachindex(parts)
        lapcof=laplacianCoeffs(parts[i])
        mu_col=collect(keys(lapcof))
        for j in eachindex(mu_col)
            if sum(mu_col[j])<N
                append!(I,i)
                append!(J,dict_ind[mu_col[j]])
                append!(V,lapcof[mu_col[j]])
            end
        end
    end
    n_len=length(parts)
    sparse(I,J,V,n_len,n_len)

end

function basisTable(N::Int, xs::AbstractVector, neg::Bool)
    P   = neg ? Ultraspherical(-0.5) : Ultraspherical(3/2)
    off = neg ? 2 : 0
    Matrix(P[xs, (1 + off):(N + off)])
end

function ApproximationPointError(N,U::Vector{Float64},u,X::Vector{Float64})
    parts=reduce(vcat, Partition_n_parts(i, 3) for i=0:N-1)
    if length(U)!=length(parts)
        error("dimension mismatch")
    end
    ac=0.0
    for i in eachindex(parts)
        ac+=U[i]* SymmetricOrthogonalPolynomials.partitionToInvariantEval(X,parts[i],true)
    end
    abs(ac-u(X...))
end

function ApproximationHeatmap(N, U::Vector{Float64}, u; xs=range(-1, 1, length=75), z=0.5)
    parts = reduce(vcat, Partition_n_parts(i, 3) for i=0:N-1)
    if length(U) != length(parts)
        error("dimension mismatch")
    end
    fs = [SymmetricOrthogonalPolynomials.partitionToInvariantFunction(parts[i], true) for i in eachindex(parts)]
    uap(X) = sum(U[i] * fs[i](X) for i in eachindex(fs))
    err = [abs(uap([x, y, z]) - u(x, y, z)) for y in xs, x in xs]
    scale = 1e7   # pick so values are O(1)
    heatmap(xs, xs, err .* scale; colorbar_title = "error (×1e-7)")
end

function errorConvergence(U::Vector{Float64}, dim::Int, N::Int; X = 2 .* rand(dim) .- 1)
    T     = basisTable(N, X, true)                        # dim × N, T[j,k+1] = basis_k(X[j]), computed ONCE
    u     = uExact(dim)
    u_X   = u(X)
    uap   = 0.0
    err   = Float64[]; degr = Int[]
    idx   = 0
    for deg in 0:N-1
        for part in Partition_n_parts(deg, dim)          # same order as assembleA!
            idx += 1
            c = U[idx]
            if c!=0.0
                for perm in multiset_permutations(part, dim) # lazy, no collect
                    p = c
                    for j in 1:dim
                        p *= T[j, perm[j] + 1]
                    end
                    uap += p
                end
            end
        end
        push!(err, abs(uap - u_X))
        push!(degr, deg)
        err[end]>2.22e-16 || break
    end
    (err, degr)
end

function ConvergencePlot(dims::Vector{Int},N::Int)
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel="Polynomial degree", ylabel="Error",
              yscale=log10, title="Approximation convergence")
    err=Float64[]
    degr=Int64[]
    degm=0
    mark=[:circle,:utriangle,:rect,:hexagon]
    aux=1
    for i in dims
        A=assembleA(N,i)
        F=assembleF(N,i)
        U=A'\F
        U[abs.(U).<1e-16].=0.0
        err,degr=errorConvergence(U,i,N)
        degm= degr[end]>degm ? degr[end] : degm
        scatterlines!(ax, degr, err; marker=mark[aux], label="dim=$i", strokewidth=1, strokecolor=:black)
        println("$i done")
        println(err[end])
        aux+=1
    end
    lines!(ax, [0,degm+1], [2.22e-16,2.22e-16]; label="Machine precission threshold")
    CairoMakie.xlims!(ax, 0, degm+1)
    axislegend(ax)
    N=N-1
    save("ConvergencePlot.png", fig)
    fig
end

function getU(N::Int, dim::Int)
    A=assembleA(N,dim)
    F=assembleF(N,dim)
    vec(A'\F)
end

function DOFforDegree(N::Int,dim::Int)
    dof=[0]
    for i in 0:N-1
        part=Partition_n_parts(i,dim)
        append!(dof,length(part)+dof[end])
    end
    (0:N-1,dof[2:end])
end

function DofvsDegPlo(dims::Vector{Int},N::Int)
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel="Polynomial degree", ylabel="Degrees of freedom",
              yscale=log10, title="Degrees of freedom vs Degree of polynomial")
    mark=[:circle,:utriangle,:rect,:hexagon]
    aux=1
    for i in dims
        (degr,dof)=DOFforDegree(N,i)
        scatterlines!(ax, degr, dof; marker=mark[aux], label="dim=$i", strokewidth=1, strokecolor=:black)
        aux+=1
    end
    CairoMakie.xlims!(ax, 0, N)
    axislegend(ax, position=:rb)
    save("DofvsDeg Poisson.png", fig)
    fig
end
ConvergencePlot([3,5,8,10],30)
DofvsDegPlo([3,5,8,10],30)