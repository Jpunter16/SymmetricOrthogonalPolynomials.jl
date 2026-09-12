using ClassicalOrthogonalPolynomials, Combinatorics, SparseArrays, LinearAlgebra, CairoMakie
using SymmetricOrthogonalPolynomials, Printf

function laplacianCoeffs(α::Vector{Int})
    SymmetricOrthogonalPolynomials.laplacianCoeffs(α, false)
end

# 3-term connection -W_n = ...
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

function assembleOp(N::Int, dim::Int, coeffFun)
    parts = reduce(vcat, Partition_n_parts(i, dim) for i = 0:(N-1))
    dict_ind = makeDictInd(parts)
    I = Int[]
    J = Int[]
    V = Float64[]
    for i in eachindex(parts)
        for (μ, w) in coeffFun(parts[i])
            sum(μ) < N || continue
            push!(I, i)
            push!(J, dict_ind[μ])
            push!(V, Float64(w))
        end
    end
    n = length(parts)
    sparse(I, J, V, n, n), parts
end

# ---------- initial condition: u0 = ∏ sinpi(x_i), rank-1 in the W basis ----------
function initialCoeffs(N::Int, dim::Int, parts)
    P = Ultraspherical(-0.5)
    xg = ClassicalOrthogonalPolynomials.grid(Ultraspherical(3/2)[:, 1:N])
    Vw = P[xg, 3:(N+2)]                       # columns W_0 … W_{N-1}
    c = Vw \ sinpi.(xg)                    # 1D coefficients of sinpi in W
    [prod(c[k+1] for k in p) for p in parts]   # rank-1, orbit factor cancels
end

# ---------- Crank–Nicolson march ----------
function heatSolve(N::Int, dim::Int, T::Float64, nsteps::Int)
    A, parts = assembleOp(N, dim, laplacianCoeffs)
    S, _ = assembleOp(N, dim, conversionCoeffs)
    dt = T / nsteps
    LHS = lu(S' - (dt/2) * A')
    RHS = S' + (dt/2) * A'
    u = initialCoeffs(N, dim, parts)
    for _ = 1:nsteps
        u = LHS \ (RHS * u)
    end
    u, parts
end

# ---------- error at final time (fast: 1D basis table, no quasi-indexing in loops) ----------
function finalError(u, parts, dim, T; xs=2 .* rand(1) .- 1)
    Tb = Matrix(Ultraspherical(-0.5)[xs, 3:(maximum(maximum, parts)+3)])
    grid = collect(Iterators.product(ntuple(_ -> 1:length(xs), dim)...))
    e = 0.0
    for I in grid
        ap = 0.0
        for (i, part) in enumerate(parts)
            for perm in multiset_permutations(part, dim)
                p = u[i]
                for j = 1:dim
                    p *= Tb[I[j], perm[j]+1]
                end
                ap += p
            end
        end
        ex = prod(sinpi(xs[I[j]]) for j = 1:dim) * exp(-dim * pi^2 * T)
        e = max(e, abs(ap - ex))
    end
    e
end

#= ---------- demo: O(dt²) convergence ----------
N, dim, T = 16, 3, 0.05
err=Float64[]
@time for m in (10, 20, 40, 80,160)
    u, parts = heatSolve(N, dim, T, m)
    append!(err, finalError(u, parts, dim, T))
end

for i in eachindex(err)
    if i+1<=length(err)
    println("Convergence: $(log2(err[i]/err[i+1]))")
    end
end

=#

function errorConvergence(dim::Int, N::Int, T::Float64, nsteps::Int)
    err = Float64[]
    degr = Int[]
    xp=2 .* rand(1) .- 1
    for deg = 1:N
        u, parts=heatSolve(deg, dim, T, nsteps)
        push!(err, finalError(u, parts, dim, T; xs=xp))
        push!(degr, deg-1)
        err[end]>2.22e-16 || break
    end
    (err, degr)
end

function ConvergencePlot(dims::Vector{Int}, N::Int, T::Float64, nsteps::Int)
    fig = Figure()
    h=@sprintf("%.2e", T/nsteps)
    ax = Axis(
        fig[1, 1];
        xlabel="Polynomial degree",
        ylabel="Error",
        yscale=log10,
        title="Approximation convergence, h=" * h * ",T=$T",
    )
    err=Float64[]
    degr=Int64[]
    degm=0
    mark=[:circle, :utriangle, :rect, :hexagon]
    aux=1
    for i in dims
        (err, degr)=errorConvergence(i, N, T, nsteps)
        degm = degr[end]>degm ? degr[end] : degm
        scatterlines!(
            ax,
            degr,
            err;
            marker=mark[aux],
            label="dim=$i",
            strokewidth=1,
            strokecolor=:black,
        )
        println("$i done")
        println(err[end])
        aux += 1
    end
    #lines!(ax, [0,degm+1], [2.22e-16,2.22e-16]; label="Machine precission")
    CairoMakie.xlims!(ax, 0, degm+1)
    Legend(fig[1, 2], ax)
    N=N-1
    save("ConvergencePlotHeat.png", fig)
    fig
end


ConvergencePlot([3, 5, 8, 10], 24, 0.05, 80)
