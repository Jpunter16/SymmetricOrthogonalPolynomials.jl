
using ClassicalOrthogonalPolynomials, Combinatorics, SymmetricOrthogonalPolynomials, SparseArrays, Plots, Printf


#f(x,y,z)=2*(y^2*z^2-z^2-y^2+1+  x^2*z^2-z^2-x^2+1  +x^2*y^2-x^2-y^2+1)

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

function assembleF(N::Int,dim::Int)
    f=RHS(dim)
    Q = Ultraspherical(3/2)
    xg = ClassicalOrthogonalPolynomials.grid(Q[:, 1:N])                # N transform points
    Vm = Q[xg, 1:N]                     # N×N generalized Vandermonde
    C = [f(collect(x)) for x in Iterators.product(ntuple(_ -> xg, dim)...)]
    for i in 1:dim
        C = mapslices(v -> Vm \ v, C; dims = i)
    end  

    fhat = Dict{Vector{Int}, Float64}()

    for i in 0:N-1
        for μ in Partition_n_parts(i, dim)      # your Partition_n_parts enumerator
            orb = multiset_permutations(μ, dim)
            fhat[μ] = sum(C[(p .+1)...] for p in orb) / SymmetricOrthogonalPolynomials.orbitSize(μ)
        end
    end
    assembleF(N,fhat,dim)
end
#=
f(x,y,z)=-3* pi^2 *sinpi(x)*sinpi(y)*sinpi(z)

Q = Ultraspherical(3/2)
N = 15                                  # tensor degrees 0:N-1 per variable

# --- tensor transform: values on grid -> C^{(3/2)} tensor coefficients ---
xg = ClassicalOrthogonalPolynomials.grid(Q[:, 1:N])                # N transform points
Vm = Q[xg, 1:N]                     # N×N generalized Vandermonde
C  = [f(xi,xj,xk) for xi in xg, xj in xg, xk in xg]
for dim in 1:3
    global C = mapslices(v -> Vm \ v, C; dims = dim)
end                                  # C[i,j,k] = c_{(i-1,j-1,k-1)}

# --- read off invariant coefficients at sorted indices ---
fhat = Dict{Vector{Int}, Float64}()


for i in 0:N-1

for μ in Partition_n_parts(i, 3)      # your Partition_n_parts enumerator
    orb = multiset_permutations(μ, 3)
    fhat[μ] = sum(C[p[1]+1, p[2]+1, p[3]+1] for p in orb) / SymmetricOrthogonalPolynomials.orbitSize(μ)
end


end

=#

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

function assembleF(N,fhat::Dict{Vector{Int}, Float64}, dim)
    parts=reduce(vcat, Partition_n_parts(i, dim) for i=0:N-1)
    dict_ind=makeDictInd(parts)
    F=zeros((length(parts),1))
    for i in keys(fhat)
        F[dict_ind[i],1]=fhat[i]
    end
    F
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

#ApproximationPointError(N,U,u,[0.1,0.2,0.3])

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

function MaxErrorGrid(uap, u_grid::Array{Float64}, grid_pts)
    uap_grid=[uap(x) for x in grid_pts]
    maximum(abs,(uap_grid.-u_grid))
end

function errorConvergence(U::Vector{Float64},dim::Int, N::Int; xg=range(-1,1,length=20))
    parts = Vector{Vector{Int64}}()
    err=Float64[]
    u=uExact(dim)
    fs=[]
    Dof=[0]

    # Precompute the evaluation grid and the exact-solution values
    grid_pts = [collect(x) for x in Iterators.product(ntuple(_ -> xg, dim)...)]
    u_grid = [u(x) for x in grid_pts]

    # Define uapprox ONCE: it closes over `fs` and `U` by reference, so as
    # `fs` grows via append! below, calling uapprox automatically sees the
    # new terms — no need to redefine the closure (and pay recompilation
    # cost) on every iteration.
    uapprox(X) = sum(U[idx]*fn(X) for (idx, fn) in enumerate(fs))

    for i in 1:N
        new_parts=Partition_n_parts(i-1,dim)
        append!(parts, new_parts)
        append!(fs, [partitionToInvariantFunction(part,true) for part in new_parts])

        # One error evaluation per degree
        append!(err, MaxErrorGrid(uapprox,u_grid,grid_pts))
        err
        append!(Dof,length(parts))
    end
    (err,Dof[2:end])
end

function ConvergencePlot(dims::Vector{Int},N::Int)
    p=plot()
    local err, Dof
    for i in dims
        A=assembleA(N,i)
        F=assembleF(N,i)
        U=vec(A'\F)
        err,Dof=errorConvergence(U,i,N)
        plot!(p,Dof,err)
        err[end]>2.22e-16 || break
    end
    print(err[end])
    p
end



function ConvergencePlotToMachinePrecision(dims::Vector{Int})
    p=plot()
    N=1
    for i in dims
        while true
            A=assembleA(N,i)
            F=assembleF(N,i)
            U=vec(A'\F)
            err,Dof=errorConvergence(U,i,N)
            err[end]>2.22e-16 || break
            N+=1
        end
        plot!(p,Dof,err)
    end 
    p
end

ConvergencePlot([3,4,5],25)

#ConvergencePlotToMachinePrecision([3])
