function computeDeltaComb(dim::Int)
    delta=[-2,0,2]
    delta_comb=vec(collect(Iterators.product(fill(delta,dim)...)))
    filter!(t->any(s->s==0,t),delta_comb)
    delta_comb
end

function nonZeroIndex(p,n_sum,delta_comb)
    p_adyacency=[p.+ d for d in delta_comb ]
    p_ady_vec=unique(p_adyacency)
    filter!(t->all(s->s>=0,t),p_ady_vec)
    filter!(t->sum(t)<=n_sum,p_ady_vec)
end

function orderPartitions(v)
    Part_dict=Dict{Vector{Int},Vector{Vector{Int}}}()
    for x in v
        key=sort(x, rev=true)
        if haskey(Part_dict, key)
            push!(Part_dict[key], x)
        else
            Part_dict[key]=[x]
        end
    end
    Part_dict
end

function getAdjacentPartitions(p,n,delta_comb)
    perm_uniq=collect(multiset_permutations(p,length(p)))
    dict=Dict{Vector{Int},Dict{Vector{Int},Vector{Vector{Int}}}}()
    for i in eachindex(perm_uniq)
        dict[perm_uniq[i]]=orderPartitions(nonZeroIndex(perm_uniq[i],n,delta_comb))
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
    aux=0.0
    #stiff_mat=(P'*diff(P,2))
    #mass_mat=(P'*P)
    for i in eachindex(same_e)
        aux_i=1.0 
        for j in eachindex(inner_prod_ind)
            if j== same_e[i]
                aux_i*=stiff_mat[inner_prod_ind[j]...]
            else
                aux_i*=mass_mat[inner_prod_ind[j]...]
            end
        end
        #if aux_i==0
        #   error("calculated 0 entry")
        #end
        aux+=aux_i
    end
    aux
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
    I_idx=Int[]
    J_idx=Int[]
    V=Float64[]
    dic=makeDictInd(partitions)
    for i=1:n_len
        dict_part=getAdjacentPartitions(partitions[i],n,delta_comb)
        adjacent_part=DistinctDictKeys(dict_part)
        for j in eachindex(adjacent_part)
            aux=0.0
            for p_perm in keys(dict_part)
                if haskey(dict_part[p_perm],adjacent_part[j])
                    for j_perm in dict_part[p_perm][adjacent_part[j]]
                        aux+=StiffMassSparseMult(p_perm,j_perm,mass_mat,stiff_mat)
                    end
                end
            end
            j_ind=dic[adjacent_part[j]]
            N_dif_p=length(unique(partitions[i]))
            N_dif_q=length(unique(adjacent_part[j]))
            push!(I_idx,i)
            push!(J_idx,j_ind)
            push!(V,aux*(1/sqrt(N_dif_p*N_dif_q)))
        end
    end
    sparse(I_idx,J_idx,V,n_len,n_len)
end
