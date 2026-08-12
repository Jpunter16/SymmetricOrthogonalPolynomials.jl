using ClassicalOrthogonalPolynomials, Combinatorics

function computeDeltaComb(dim::Int)
    delta=[-2,0,2]
    delta_comb=vec(collect(Iterators.product(fill(delta,dim)...)))
    filter!(t->any(s->s==0,t),delta_comb)
    delta_comb
end

function nonZeroIndeX(p,n_sum,delta_comb)
    p_adyacency=[p.+ d for d in delta_comb ]
    p_ady_vec=unique(p_adyacency)
    filter!(t->all(s->s>=0,t),p_ady_vec)
    filter!(t->sum(t)<=n_sum,p_ady_vec)
end

function orderPartitions(v)
    Part_dict=Dict{Vector{Int},Vector{Vector{Int}}}()
    for X in v
        key=sort(X, rev=true)
        if haskey(Part_dict, key)
            push!(Part_dict[key], X)
        else
            Part_dict[key]=[X]
        end
    end
    Part_dict
end

function getAdjacentPartitions(p,n,delta_comb)
    perm_uniq=collect(multiset_permutations(p,length(p)))
    dict=Dict{Vector{Int},Dict{Vector{Int},Vector{Vector{Int}}}}()
    for i in eachindeX(perm_uniq)
        dict[perm_uniq[i]]=orderPartitions(nonZeroIndeX(perm_uniq[i],n,delta_comb))
    end
    dict
end

function DistinctDictKeys(dict)
    keys_set=Set{Vector{Int}}()
    for key in keys(dict)
        union!(keys_set,keys(dict[key]))
    end
    collect(keys_set)
end

function StiffMassSparseMult(p,q, mass_mat,stiff_mat)
    same_e=findall(p .==q)
    inner_prod_ind=collect(zip(p.+1,q.+1))
    auX=0.0
    #stiff_mat=(P'*diff(P,2))
    #mass_mat=(P'*P)
    for i in eachindeX(same_e)
        auX_i=1.0 
        for j in eachindeX(inner_prod_ind)
            if j== same_e[i]
                auX_i*=stiff_mat[inner_prod_ind[j]...]
            else
                auX_i*=mass_mat[inner_prod_ind[j]...]
            end
        end
        #if auX_i==0
        #   error("calculated 0 entry")
        #end
        auX+=auX_i
    end
    auX
end

function makeDictInd(part_v::Vector{Vector{Int64}})
    dic=Dict{Vector{Int},Int}()
    for (i,el) in enumerate(part_v)
        dic[el]=i
    end
    dic
end

function getLaplacianClosedForm(P,n::Int, dim::Int)
    partitions=[part for i=0:n for part in Partition_n_parts(i,dim)]
    n_len=length(partitions)
    stiff_full = P' * diff(P, 2)
    stiff_mat = Diagonal([stiff_full[i,i] for i in 1:n_len])
    mass_mat = (P' * P)[1:n_len,1:n_len]
    delta_comb = computeDeltaComb(dim)
    I_idX=Int[]
    J_idX=Int[]
    V=Float64[]
    dic=makeDictInd(partitions)
    for i=1:n_len
        dict_part=getAdjacentPartitions(partitions[i],n,delta_comb)
        adjacent_part=DistinctDictKeys(dict_part)
        for j in eachindeX(adjacent_part)
            auX=0.0
            for p_perm in keys(dict_part)
                if haskey(dict_part[p_perm],adjacent_part[j])
                    for j_perm in dict_part[p_perm][adjacent_part[j]]
                        auX+=StiffMassSparseMult(p_perm,j_perm,mass_mat,stiff_mat)
                    end
                end
            end
            j_ind=dic[adjacent_part[j]]
            N_dif_p=length(unique(partitions[i]))
            N_dif_q=length(unique(adjacent_part[j]))
            push!(I_idX,i)
            push!(J_idX,j_ind)
            push!(V,auX*(1/sqrt(N_dif_p*N_dif_q)))
        end
    end
    sparse(I_idX,J_idX,V,n_len,n_len)
end

### Other Way

# 3-term connection -W_n = ...
splus1(n) = -1 // ((2n + 3) * (2n + 5))
s0(n) =   2 // ((2n + 1) * (2n + 5))
sminus1(n) = -1 // ((2n + 1) * (2n + 3))

# |Stab(μ)| = ∏ mult! (zeros count as slots)
stab(μ) = prod(factorial(count(==(v), μ)) for v in unique(μ))

# admissible (up,stay,down), up+stay+down = m: shift 2/0/-2 needs degree >= 0/0/2
comps(v, m) = [(up, stay, m - up - stay) for up in 0:m for stay in 0:(m - up)
               if (m - up - stay == 0 || v >= 2)]


orbitSize(α::Vector{Int}) =factorial(length(α)) ÷ stab(α)

function laplacianCoeffs(α::Vector{Int}, normalized::Bool = false)
    vals = sort(unique(α); rev = true)                 # distinct parts v₁ > … > v_r
    mult = Dict(v => count(==(v), α) for v in vals)    # multiplicities m_u

    buckets = Dict{Vector{Int}, Float64}()

    for (t, vt) in enumerate(vals)
        lists = [comps(v, mult[v] - (u == t)) for (u, v) in enumerate(vals)]
        for combo in Iterators.product(lists...)
            w = 1 / 1
            landing = Int[vt ]
            for (u, (up, stay, down)) in enumerate(combo)
                v = vals[u]
                w *= splus1(v)^up * s0(v)^stay * sminus1(v)^down /
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
        μ_norm= Dict(μ => 1/sqrt(orbitSize(μ)) for μ in collect(keys(buckets)) )
        Dict(μ => α_norm* μ_norm[μ]* stab(μ) * w for (μ, w) in buckets if !iszero(w))
    else
        Dict(μ => -stab(μ) * w for (μ, w) in buckets if !iszero(w))
    end
end


#laplacianCoeffs([9,8,7,6,5,4,3,2,0])

function partitionToInvariantEval(X::Vector{Float64},α::Vector{Int}, neg::Bool)
    n=length(α)
    if neg
    P=Ultraspherical(-0.5)[:, 3:end]
    else
       P=Ultraspherical(1.5) 
    end
    perms = collect(multiset_permutations(α, n))
    auX=0.0
    for i in eachindex(perms)
        auXm=1.0
        for j in 1:n
            auXm*= P[X[j], perms[i][j] + 1]
        end
        auX += auXm
    end
    auX
end

function partitionToInvariantFunction(α::Vector{Int}, neg::Bool)
    n=length(α)
    if neg
    P=Ultraspherical(-0.5)[:, 3:end]
    else
       P=Ultraspherical(1.5) 
    end
    perms = collect(multiset_permutations(α, n))
    return X::AbstractVector -> begin
        @assert length(X) == n "Expected vector of length $n"
        auX=0.0
        for i in eachindex(perms)
            auXm=1.0
            for j in 1:n
                auXm*= P[X[j], perms[i][j] + 1]
            end
            auX += auXm
        end
        auX
    end
end

function laplacianFiniteDiff3D(X::Vector{Float64},α::Vector{Int},h::Float64)
    u(x)=partitionToInvariantEval(x,α,true)
    (u([X[1]+h,X[2],X[3]])+u([X[1]-h,X[2],X[3]])+u([X[1],X[2]+h,X[3]])
    +u([X[1],X[2]-h,X[3]])+u([X[1],X[2],X[3]+h])+u([X[1],X[2],X[3]-h])
    -6*u([X[1],X[2],X[3]]))/(h^2)
end

function LapCoefftoEval(X::Vector{Float64},α::Vector{Int})
    dic=laplacianCoeffs(α,false)
    aux=0.0
    for i in collect(keys(dic))
        aux+=dic[i]* partitionToInvariantEval(X,i,false)
    end
    aux
end

X=[0.1,0.2,0.3]
alpha=[3,2,1]
h=0.001
laplacianFiniteDiff3D(X,alpha,h)-LapCoefftoEval(X,alpha)