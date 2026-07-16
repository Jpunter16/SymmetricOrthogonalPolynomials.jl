function nonZeroIndex(p,n_sum)
    n=length(p)
    delta=[-2,0,2]
    delta_comb=vec(collect(Iterators.product(fill(delta,n)...)))
    filter!(t->any(s->s==0,t),delta_comb)
    p_adyacency=[p.+ collect(d) for d in delta_comb ]
    p_ady_vec=unique(vec(p_adyacency))
    filter!(t->all(s->s>=0,t),p_ady_vec)
    filter!(t->sum(t)<=n_sum,p_ady_vec)
end

function orderPartitions(v)
    b=copy(v)
    taken=falses(length(b))
    Part_dict=Dict()
    for i=eachindex(v)
        if !taken[i]
            perm=sort(b[i], rev=true)
            reordered_perm=[b[i]]
            for j=(i+1):length(b)
                if sort(b[j],rev=true)==perm
                    taken[j]=true
                    push!(reordered_perm, b[j])
                end
            end
            Part_dict[perm]=reordered_perm
        end
    end
    Part_dict
end

function addDictionary(dict_o::Dict,dict_a::Dict)
    dict_or=deepcopy(dict_o)
    dict_add=deepcopy(dict_a)
    for i in keys(dict_add) 
        if haskey(dict_or, i)
            push!(dict_or[i],dict_add[i]...)
        else
            dict_or[i]=dict_add[i]
        end
    end
    dict_or
end

function getAdjacentPartitions(p,n)
    perm=collect(permutations(p))
    perm_uniq=unique(perm)
    dict=Dict()
    for i in eachindex(perm_uniq)
        dict[perm_uniq[i]]=orderPartitions(nonZeroIndex(perm_uniq[i],n))
    end
    dict
end

function DistinctDictKeys(dict)
    unique([keys_in for key in keys(dict) for keys_in in keys(dict[key])])
end

function StiffMassSparseMult(p,q,P)
    same_e=findall(p .==q)
    inner_prod_ind=collect(zip(p.+1,q.+1))
    aux=0.0
    stiff_mat=(P'*diff(P,2))
    mass_mat=(P'*P)
    for i in eachindex(same_e)
        aux_i=1.0 
        for j in eachindex(inner_prod_ind)
            if j== same_e[i]
                aux_i*=stiff_mat[inner_prod_ind[j]...]
            else
                aux_i*=mass_mat[inner_prod_ind[j]...]
            end
        end
        aux+=aux_i
    end
    aux
end

function getLaplacianClosedForm(P,n::Int, dim::Int)
    partitions=[part for i=0:n for part in Partition_n_parts(i,dim)]
    n_len=length(partitions)
    Δ=spzeros(n_len,n_len)
    for i=1:n_len
        dict_part=getAdjacentPartitions(partitions[i],n)
        adjacent_part=DistinctDictKeys(dict_part)
        for j in eachindex(adjacent_part)
            aux=0.0
            for p_perm in keys(dict_part)
                if haskey(dict_part[p_perm],adjacent_part[j])
                    aux+=StiffMassSparseMult(p_perm,adjacent_part[j],P)
                end
            end
            #print(adjacent_part[j])
            Δ[i,findfirst(x->x==adjacent_part[j], partitions)]=aux
        end
    end
    Δ
end
