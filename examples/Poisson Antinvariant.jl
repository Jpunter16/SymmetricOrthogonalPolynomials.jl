using LinearAlgebra, ClassicalOrthogonalPolynomials, SymmetricOrthogonalPolynomials, CairoMakie, SparseArrays, Combinatorics

uExact(xs)  = det([sinpi(j*xi) for xi in xs, j in 1:length(xs)])

function assembleF(N, fhat::Dict{Vector{Int},Float64}, dim)
    parts=reduce(vcat, Partition_diff_n_parts(i, dim) for i=0:(N-1))
    dict_ind=makeDictInd(parts)
    F=zeros((length(parts)))
    for i in keys(fhat)
        F[dict_ind[i]]=fhat[i]
    end
    F
end

function assembleF_det(N::Int, dim::Int)
    N=Int((dim*(dim+1))/2)+N

    S = dim*(dim+1)*(2dim+1)/6
    scale = -S*pi^2                       # f = Δu = scale * u

    C = Ultraspherical(3/2)
    x = axes(C, 1)
    cs = [C \ sinpi.(j .* x) for j in 1:dim]   # coefficients of sin(π j x)

    # coefficient lookup, 0 beyond the computed length
    coef(j, k) = k + 1 <= length(cs[j]) ? cs[j][k+1] : 0.0

    fhat = Dict{Vector{Int},Float64}()
    for i in 0:(N-1)
        for μ in Partition_diff_n_parts(i, dim)
            allunique(μ) || continue      # repeated index => det = 0
            A = [coef(j, k) for k in μ, j in 1:dim]   # row ↔ μ_i, column ↔ j
            fhat[μ] = scale * det(A)
        end
    end
    assembleF(N, fhat, dim)               # your existing Dict method
end

function assembleA(N::Int, dim::Int)
    N=Int((dim*(dim+1))/2)+N
    parts = reduce(vcat, Partition_diff_n_parts(i, dim) for i=0:(N-1))
    dict_ind = makeDictInd(parts)
    I = Int[]
    J = Int[]
    V = Float64[]
    for i in eachindex(parts)
        lapcof = laplacianCoeffsAlter(parts[i])
        mu_col = collect(keys(lapcof))
        for j in eachindex(mu_col)
            if sum(mu_col[j])<N
                append!(I, i)
                append!(J, dict_ind[mu_col[j]])
                append!(V, lapcof[mu_col[j]])
            end
        end
    end
    n_len = length(parts)
    sparse(I, J, V, n_len, n_len)
end

function basisTable(N::Int, xs::AbstractVector, neg::Bool)
    P = neg ? Ultraspherical(-0.5) : Ultraspherical(3/2)
    off = neg ? 2 : 0
    Matrix(P[xs, (1+off):(N+off)])
end

function errorConvergence(U::AbstractVector{<:Real}, dim::Int, N::Int; X=2 .* rand(dim) .- 1)
    N=Int((dim*(dim+1))/2)+N
    T = basisTable(N, X, true)
    u_X = uExact(X)
    uap = 0.0
    err = Float64[]
    degr = Int[]
    idx = 0
    for deg in 0:(N-1)
        for part in Partition_diff_n_parts(deg, dim)
            idx += 1
            c = U[idx]
            c == 0.0 && continue
            uap += c * det([T[j, k+1] for j in 1:dim, k in part])
        end
        if deg>= dim*(dim+1)/2
            push!(err, abs(uap - u_X))
            push!(degr, deg-Int((dim*(dim+1))/2))
        end
        if length(err)>=1
            err[end] > 2.22e-16 || break
        end
    end
    (err, degr)
end

function ConvergencePlot(dims::Vector{Int}, N::Vector{Int})
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Polynomial degree", ylabel="Error",
        yscale=log10, title="Approximation convergence")
    err=Float64[]
    degr=Int64[]
    degm = 0
    mark = [:circle, :utriangle, :rect,:pentagon, :hexagon]
    aux = 1
    for i in eachindex(dims)
        A = assembleA(N[i], dims[i])
        F = assembleF_det(N[i], dims[i])
        U = A'\F
        err, degr = errorConvergence(U, dims[i], N[i])
        degm = degr[end] > degm ? degr[end] : degm
        scatterlines!(ax, degr, err; marker=mark[aux], label="dim=$(dims[i])", strokewidth=1, strokecolor=:black)
        println("$(dims[i]) done")
        println(err[end])
        aux+=1
    end
    lines!(ax, [0, degm+1], [2.22e-16, 2.22e-16]; label="Machine precission threshold")
    CairoMakie.xlims!(ax, 0, degm+1)
    axislegend(ax, position=:lb)
    save("ConvergencePlot.png", fig)
    fig
end

function getU(N::Int, dim::Int)
    A=assembleA(N, dim)
    F=assembleF_det(N, dim)
    vec(A'\F)
end

function DOFforDegree(N::Int, dim::Int)
    dof=[0]
    for i in 0:(N-1)
        part=Partition_diff_n_parts(i, dim)
        append!(dof, length(part)+dof[end])
    end
    (0:(N-1), dof[2:end])
end

function DofvsDegPlo(dims::Vector{Int}, N::Vector{Int})
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="Polynomial degree", ylabel="Degrees of freedom",
        yscale=log10, title="Degrees of freedom vs Degree of polynomial")
    mark=[:circle, :utriangle, :rect, :hexagon]
    aux=1
    for i in eachindex(dims)
        (degr, dof)=DOFforDegree(N[i], dims[i])
        scatterlines!(ax, degr, dof; marker=mark[aux], label="dim=$(dims[i])", strokewidth=1, strokecolor=:black)
        aux+=1
    end
    CairoMakie.xlims!(ax, 0, N[i])
    axislegend(ax, position=:rb)
    save("DofvsDeg Poisson AntiInvariant.png", fig)
    fig
end

dims=[3,4, 5,6]
ConvergencePlot(dims,30*ones(Int,length(dims)))