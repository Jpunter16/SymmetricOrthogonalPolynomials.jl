
# make a tuple corresponding to lexigraphical order
function lextuples(n)
    ret = NTuple{3,Int}[]
    for k = 1:n, j=1:k
        push!(ret, (n-k+1,k-j+1,j))
    end
    ret
end

function cubegen1(n)
    tup = lextuples(n)
    rev = map(reverse, tup)
    p = sortperm(rev; rev=true)
    N = length(tup)
    ret = zeros(Int,N,N)
    for k = 1:N
        ret[k,p[k]] = (-1)^(rev[k][2]+1)
    end
    ret
end


function cubegen2(n)
    tup = lextuples(n)
    rev = map(((a,b,c),) -> (a,c,b), tup)
    p = sortperm(rev; rev=true)
    N = length(tup)
    ret = zeros(Int,N,N)
    s = (-1)^(n+1)
    for k = 1:N
        ret[k,p[k]] = s
    end
    ret
end

function cubegen3(n)
    tup = lextuples(n)
    rev = map(reverse, tup)
    p = sortperm(rev; rev=true)
    N = length(tup)
    ret = zeros(Int,N,N)
    s = (-1)^(n+1)
    for k = 1:N
        ret[k,p[k]] = s
    end
    ret
end

cuberepresentation(n) = Representation([cubegen1(n), cubegen2(n), cubegen3(n)])

function _cube_filter!(ret, ρ, p, n, j)
    λ = multiplicities(ρ)
    kys = sort!(collect(keys(λ)))
    ind = 0
    for k in kys
        m = hooklength(k)
        if k == p
            @assert 1 ≤ j ≤ m
            ret[StepRangeLen(ind+j, m, λ[k])] .= true
            return ret
        else
            ind += λ[k]*m
        end
    end
    ret
end

function _cube_filter!(ret, ρ, p, n)
    λ = multiplicities(ρ)
    kys = sort!(collect(keys(λ)))
    ind = 0
    for k in kys
        m = hooklength(k)
        if k == p
            ret[ind+1:ind+λ[k]*m] .= true
            return ret
        else
            ind += λ[k]*m
        end
    end
    ret
end


function cube_filter((p, s), N, j...)
    ret = zeros(Bool, binomial(N+2, N-1))
    ind = 0
    for n = 1:N
        M = sum(1:n)
        if s == isodd(n)
            _cube_filter!(view(ret, ind+1:ind+M), cuberepresentation(n), p, n, j...)
        end
        ind += M
    end
    ret
end