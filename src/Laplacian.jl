using ClassicalOrthogonalPolynomials, Combinatorics

function computeDeltaComb(dim::Int)
    delta=[-2, 0, 2]
    delta_comb=vec(collect(Iterators.product(fill(delta, dim)...)))
    filter!(t->any(s->s==0, t), delta_comb)
    delta_comb
end

function nonZeroIndeX(p, n_sum, delta_comb)
    p_adyacency=[p .+ d for d in delta_comb]
    p_ady_vec=unique(p_adyacency)
    filter!(t->all(s->s>=0, t), p_ady_vec)
    filter!(t->sum(t)<=n_sum, p_ady_vec)
end

function orderPartitions(v)
    Part_dict=Dict{Vector{Int},Vector{Vector{Int}}}()
    for X in v
        key=sort(X, rev = true)
        if haskey(Part_dict, key)
            push!(Part_dict[key], X)
        else
            Part_dict[key]=[X]
        end
    end
    Part_dict
end

function getAdjacentPartitions(p, n, delta_comb)
    perm_uniq=collect(multiset_permutations(p, length(p)))
    dict=Dict{Vector{Int},Dict{Vector{Int},Vector{Vector{Int}}}}()
    for i in eachindeX(perm_uniq)
        dict[perm_uniq[i]]=orderPartitions(nonZeroIndeX(perm_uniq[i], n, delta_comb))
    end
    dict
end

function DistinctDictKeys(dict)
    keys_set=Set{Vector{Int}}()
    for key in keys(dict)
        union!(keys_set, keys(dict[key]))
    end
    collect(keys_set)
end

function StiffMassSparseMult(p, q, mass_mat, stiff_mat)
    same_e=findall(p .== q)
    inner_prod_ind=collect(zip(p .+ 1, q .+ 1))
    auX=0.0
    for i in eachindeX(same_e)
        auX_i=1.0
        for j in eachindeX(inner_prod_ind)
            if j == same_e[i]
                auX_i*=stiff_mat[inner_prod_ind[j]...]
            else
                auX_i*=mass_mat[inner_prod_ind[j]...]
            end
        end
        auX+=auX_i
    end
    auX
end

function makeDictInd(part_v::Vector{Vector{Int64}})
    dic=Dict{Vector{Int},Int}()
    for (i, el) in enumerate(part_v)
        dic[el]=i
    end
    dic
end

#Deprecated way, see function laplacianCoeffs
function getLaplacianClosedForm(P, n::Int, dim::Int)
    partitions=[part for i = 0:n for part in Partition_n_parts(i, dim)]
    n_len=length(partitions)
    stiff_full = P' * diff(P, 2)
    stiff_mat = Diagonal([stiff_full[i, i] for i = 1:n_len])
    mass_mat = (P'*P)[1:n_len, 1:n_len]
    delta_comb = computeDeltaComb(dim)
    I_idX=Int[]
    J_idX=Int[]
    V=Float64[]
    dic=makeDictInd(partitions)
    for i = 1:n_len
        dict_part=getAdjacentPartitions(partitions[i], n, delta_comb)
        adjacent_part=DistinctDictKeys(dict_part)
        for j in eachindeX(adjacent_part)
            auX=0.0
            for p_perm in keys(dict_part)
                if haskey(dict_part[p_perm], adjacent_part[j])
                    for j_perm in dict_part[p_perm][adjacent_part[j]]
                        auX+=StiffMassSparseMult(p_perm, j_perm, mass_mat, stiff_mat)
                    end
                end
            end
            j_ind=dic[adjacent_part[j]]
            N_dif_p=length(unique(partitions[i]))
            N_dif_q=length(unique(adjacent_part[j]))
            push!(I_idX, i)
            push!(J_idX, j_ind)
            push!(V, auX*(1/sqrt(N_dif_p*N_dif_q)))
        end
    end
    sparse(I_idX, J_idX, V, n_len, n_len)
end

### Other Way

# 3-term connection formulas
splus1(n::Int) = -1 / ((2n + 3) * (2n + 5))
s0(n::Int) = 2 / ((2n + 1) * (2n + 5))
sminus1(n::Int) = -1 / ((2n + 1) * (2n + 3))


function s(option::Int, value::Int)
    if option==1
        return splus1(value)
    elseif option==0
        return s0(value)
    elseif option==-1
        return sminus1(value)
    end
end

stab(μ) = prod(factorial(count(==(v), μ)) for v in unique(μ))

# admissible (up,stay,down)
comps(v, m) = [
    (up, stay, m - up - stay) for up = 0:m for
    stay = 0:(m-up) if (m - up - stay == 0 || v >= 2)
]


orbitSize(α::Vector{Int}) = factorial(length(α)) ÷ stab(α)

function laplacianCoeffs(α::Vector{Int}, normalized::Bool = false)
    vals = sort(unique(α); rev = true)
    mult = Dict(v => count(==(v), α) for v in vals)

    buckets = Dict{Vector{Int},Float64}()

    for (t, vt) in enumerate(vals)
        lists = [comps(v, mult[v] - (u == t)) for (u, v) in enumerate(vals)]
        for combo in Iterators.product(lists...)
            w = 1 / 1
            landing = Int[vt]
            for (u, (up, stay, down)) in enumerate(combo)
                v = vals[u]
                w *=
                    splus1(v)^up * s0(v)^stay * sminus1(v)^down /
                    (factorial(up) * factorial(stay) * factorial(down))
                append!(landing, fill(v+2, up))
                append!(landing, fill(v, stay))
                append!(landing, fill(v - 2, down))
            end
            μ = sort(landing; rev = true)
            buckets[μ] = get(buckets, μ, 0 / 1) + w
        end
    end

    if normalized
        α_norm=1/sqrt(orbitSize(α))
        μ_norm = Dict(μ => 1/sqrt(orbitSize(μ)) for μ in collect(keys(buckets)))
        Dict(μ => -α_norm * μ_norm[μ] * stab(μ) * w for (μ, w) in buckets if !iszero(w))
    else
        Dict(μ => -stab(μ) * w for (μ, w) in buckets if !iszero(w))
    end
end

function laplacianCoeffsAlter(α::Vector{Int}, normalized::Bool = false)
    if length(unique(α))==length(α)
        buckets = Dict{Vector{Int},Float64}()
        d=length(α)
        for i in eachindex(α)
            for k in Iterators.product(ntuple(_ -> -1:1, d-1)...)
                k=(k[1:i-1]..., 1000, k[i:end]... ) #the 1000 value is not used in the implementation, it just corresponds to the differentiated slot that is already in the positive ultrasphericals
                mu = copy(α)
                for ki in eachindex(k)
                    if i !=ki
                        mu[ki]+=2*k[ki]
                    end
                end
                if all(t-> t>=0, mu)
                    if length(mu)==length(unique(mu)) #check valid alternating partition (all parts distinct)
                        pmu = sortperm(mu, rev=true)
                        T=0.0
                        q = invperm(sortperm(pmu))  
                        sig=levicivita(q)
                        aux=1.0
                        for ki in eachindex(k)
                            if i !=ki
                                aux*=s(k[ki],α[ki])
                            end
                        end
                        T+=sig*aux
                        buckets[mu[pmu]] = get(buckets, mu[pmu], 0 / 1) + T
                    end
                end
            end
        end
    end
    for key in keys(buckets) #change sign according to the differentiation rule
        buckets[key]=-buckets[key]
    end
    buckets
end

function TensorEval(α::Vector{Int},X::Vector{Float64}, neg::Bool)
    if neg
        C=Ultraspherical(-0.5)[:,3:end]
    else
        C=Ultraspherical(1.5)
    end
    aux=1.0
    for i in eachindex(α)
        aux*=C[X[i],α[i]+1]
    end
    aux
end

function AlterEval(α::Vector{Int},X::Vector{Float64}, neg::Bool)
    perms=collect(multiset_permutations(α))
    res=0.0
    for p in perms
        p_sort=invperm(sortperm(p, rev=true))
        res+=levicivita(p_sort)*TensorEval(p, X, neg)
    end
    res
end

function partitionToInvariantEval(X::Vector{Float64}, α::Vector{Int}, neg::Bool)
    n=length(α)
    if neg
        P=Ultraspherical(-0.5)[:, 3:end]
    else
        P=Ultraspherical(1.5)
    end
    vals = sort(unique(α); rev = true)
    mult = [count(==(v), α) for v in vals]
    orbitSumDP(P, X, vals, mult, n)
end

function partitionToInvariantFunction(α::Vector{Int}, neg::Bool)
    n=length(α)
    if neg
        P=Ultraspherical(-0.5)[:, 3:end]
    else
        P=Ultraspherical(1.5)
    end
    vals = sort(unique(α); rev = true)
    mult = [count(==(v), α) for v in vals]
    r = length(vals)
    dims = Tuple(m + 1 for m in mult)
    D = zeros(Float64, dims)
    Dnew = similar(D)
    return X::AbstractVector -> begin
        @assert length(X) == n "Expected vector of length $n"
        fill!(D, 0.0)
        D[ntuple(_->1, r)...] = 1.0
        for j = 1:n
            fill!(Dnew, 0.0)
            for idx in CartesianIndices(D)
                v = D[idx]
                v == 0.0 && continue
                for k = 1:r
                    if idx[k] <= mult[k]
                        t = Base.setindex(Tuple(idx), idx[k] + 1, k)
                        Dnew[t...] += v * P[X[j], vals[k]+1]
                    end
                end
            end
            D, Dnew = Dnew, D
        end
        D[ntuple(k->mult[k]+1, r)...]
    end
end

function laplacianFiniteDiff3D(X::Vector{Float64}, α::Vector{Int}, h::Float64)
    u(x) = partitionToInvariantEval(x, α, true)
    (
        u([X[1]+h, X[2], X[3]]) +
        u([X[1]-h, X[2], X[3]]) +
        u([X[1], X[2]+h, X[3]]) +
        u([X[1], X[2]-h, X[3]]) +
        u([X[1], X[2], X[3]+h]) +
        u([X[1], X[2], X[3]-h]) - 6*u([X[1], X[2], X[3]])
    )/(h^2)
end

function laplacianFiniteDiff3DAlter(X::Vector{Float64}, α::Vector{Int}, h::Float64)
    u(x) = AlterEval(α, x, true)
    (
        u([X[1]+h, X[2], X[3]]) +
        u([X[1]-h, X[2], X[3]]) +
        u([X[1], X[2]+h, X[3]]) +
        u([X[1], X[2]-h, X[3]]) +
        u([X[1], X[2], X[3]+h]) +
        u([X[1], X[2], X[3]-h]) - 6*u([X[1], X[2], X[3]])
    )/(h^2)
end

function LapCoeffToAlter(α::Vector{Int64},X::Vector{Float64})
    dic=laplacianCoeffsAlter(α)
    aux=0.0
    for i in keys(dic) 
        aux+=dic[i]*AlterEval(i, X, false)
    end
    aux
end

function LaplacianErrorAlter(α::Vector{Int64},X::Vector{Float64}, h::Float64)
    abs(LapCoeffToAlter(α, X)-laplacianFiniteDiff3DAlter(X, α, h))
end

function LapCoefftoEval(X::Vector{Float64}, α::Vector{Int})
    dic=laplacianCoeffs(α, false)
    aux=0.0
    for i in collect(keys(dic))
        aux+=dic[i] * partitionToInvariantEval(X, i, false)
    end
    aux
end

function LaplacianMat(N::Int, dim::Int)
    parts=reduce(vcat, Partition_n_parts(i, dim) for i = 0:(N-1))
    dict_ind = makeDictInd(parts)
    n_len=length(parts)
    diag_11=[2/(2*n+3) for n = 0:(N-1)]
    I = Int[]
    J = Int[]
    V = Float64[]
    for i in eachindex(parts)
        lapcof=laplacianCoeffs(parts[i])
        for j in keys(lapcof)
            if sum(j)<=N-1
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

for h in [0.1,0.01,0.001]
    println(LaplacianErrorAlter([5,4,3],[0.1,0.2,0.3],h))
end