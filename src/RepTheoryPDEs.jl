using NumericalRepresentationTheory, SparseArrays, BlockArrays
export cubeperm_Q, cubeperm_inds, cubeperm_filter, cubevec2tens
import SymmetricOrthogonalPolynomials: lextuples
import NumericalRepresentationTheory:partitions

function cubeperm_gen1(n)
    tup = lextuples(n)
    rev = map(((a,b,c),) -> (b,a,c), tup)
    p = sortperm(rev; rev=true)
    N = length(tup)
    ret = zeros(Int,N,N)
    for k = 1:N
        ret[k,p[k]] = 1
    end
    ret
end
function cubeperm_gen2(n)
    tup = lextuples(n)
    rev = map(((a,b,c),) -> (a,c,b), tup)
    p = sortperm(rev; rev=true)
    N = length(tup)
    ret = zeros(Int,N,N)
    for k = 1:N
        ret[k,p[k]] = 1
    end
    ret
end

cubeperm_representation(n) = Representation([cubeperm_gen1(n), cubeperm_gen2(n)])


function cubeperm_Q(N)
    Qs = Matrix{Float64}[]
    for n = 1:N
        @show n
        @time push!(Qs, blockdiagonalize(cubeperm_representation(n))[2])
    end

    Q = blockdiag(sparse.(Qs)...)
end

# find all indices for irrep with partition p and SYT corresponding to j
function _cubeperm_filter!(ret, p, n, j)
    λ = multiplicities(cubeperm_representation(n)) # mults of all irreps
    kys = sort!(collect(keys(λ))) # sort partitions present
    ind = 0
    for k in kys # for each irrep
        m = hooklength(k) # dimension of irrep
        if k == p # the ones we care about
            @assert 1 ≤ j ≤ m
            ret[StepRangeLen(ind+j, m, λ[k])] .= true # only take every m row up to the given multiplicity
            return ret
        else
            ind += λ[k]*m
        end
    end
    ret
end

# find all indices for irrep with partition p
function _cubeperm_filter!(ret, p, n)
    λ = multiplicities(cubeperm_representation(n)) # find multiplicities for all irreps
    kys = sort!(collect(keys(λ))) # which irreps are present
    ind = 0
    for k in kys
        m = hooklength(k) # dim of irrep
        if k == p # irrep we care about
            ret[ind+1:ind+λ[k]*m] .= true # the next m indices correspond to this irrep
            return ret
        else
            ind += λ[k]*m # ignore these entries
        end
    end
    ret
end


function cubeperm_filter((p, s), N, j...)
    # vector of 0 or 1 determining which indices are in the irrep corresponding to partition p and SYT given by j
    # up to degree N-1
    ret = zeros(Bool, binomial(N+2, N-1))
    ind = 0
    for n = 1:N # degree n-1 polynomials
        M = sum(1:n) # dim of degree n-1 polys
        if s == isodd(n) # check if we match even/odd
            _cubeperm_filter!(view(ret, ind+1:ind+M), p, n, j...) # populate degree n-1 case
        end
        ind += M
    end
    ret
end


function cubeperm_inds(N)
    inds = Int[]
    for s in (true,false), p in NumericalRepresentationTheory.partitions(3)
        for j = 1:hooklength(p)
            append!(inds, findall(cubeperm_filter((p,s), N, j)))
        end
    end
    inds
end

function cubevec2tens(C, N)
    M = blocksize(C,1)
    X = zeros(N, N, N)
    for n = 1:M
        blk = C[Block(n)]
        inds = lextuples(n)
        for (val,ind) in zip(blk,inds)
            X[ind...] = val
        end
    end
    X
end


function S3AntiAndInvariantBasisPruner(n::Int)
    alltups = vcat(([part for i in 0:n for part in addDistinctPartition(i)[1]]))

    inds_o = findall(t -> isodd(sum(t)), alltups)
    inds_e = findall(t -> iseven(sum(t)), alltups)
    inds_parity=[inds_o;inds_e]

    inds=[1:length(alltups)...]
    inds_keep = findall(t -> t[1] >= 2 && t[2] >= 2 && t[3] >=2, alltups)

    inds_keep_parity=filter(in(inds_keep),inds_parity)

    Anti_mask=[part for i in 0:n for part in addDistinctPartition(i)[2]]
    inds_Anti=findall(Anti_mask)
    inds_Invar=filter(!in(inds_Anti),inds)
    inds_keep_Anti=filter(in(inds_keep),inds_Anti)
    inds_keep_Invar=filter(in(inds_keep),inds_Invar)
    inds_keep_Invariance=[inds_keep_Anti;inds_keep_Invar]

    inds_anti_odd=filter(in(inds_o),inds_keep_Anti)
    inds_anti_even=filter(in(inds_e),inds_keep_Anti)
    inds_Invar_odd=filter(in(inds_o),inds_keep_Invar)
    inds_Invar_even=filter(in(inds_e),inds_keep_Invar)
    inds_keep_Invariance_parity=[inds_anti_odd;inds_anti_even;inds_Invar_odd;inds_Invar_even]
    vcat(([inds_keep],[inds_keep_Invariance]),
    ([inds_keep_parity],[inds_keep_Invariance_parity]))
end

perm_representation(n)=Representation([S3perm_gen1(n),S3perm_gen2(n)])

# find all indices for irrep with partition p and SYT corresponding to j
function _perm_filter!(ret, p, n, j)
    λ = multiplicities(perm_representation(n)) # mults of all irreps
    kys = sort!(collect(keys(λ))) # sort partitions present
    ind = 0
    for k in kys # for each irrep
        m = hooklength(k) # dimension of irrep
        if k == p # the ones we care about
            @assert 1 ≤ j ≤ m
            ret[StepRangeLen(ind+j, m, λ[k])] .= true # only take every m row up to the given multiplicity
            return ret
        else
            ind += λ[k]*m
        end
    end
    ret
end

function perm_filter((p), N, j...)
    # vector of 0 or 1 determining which indices are in the irrep corresponding to partition p and SYT given by j
    # up to degree N-1
    ret = zeros(Bool, binomial(N+3,3))
    ind = 0
    for n = 1:N # degree n-1 polynomials
        M = binomial(n+1,2) # dim of degree n-1 polys
        _perm_filter!(view(ret, ind+1:ind+M), p, n, j...) # populate degree n-1 case
        ind += M
    end
    ret
end

function perm_inds(N)
    inds = Int[]
    for p in NumericalRepresentationTheory.partitions(3)
        for j = 1:hooklength(p)
            append!(inds, findall(perm_filter((p), N, j)))
        end
    end
    inds
end