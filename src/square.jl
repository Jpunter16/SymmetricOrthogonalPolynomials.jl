
# make a tuple corresponding to lexigraphical order
function lextuples2(n)
    ret = NTuple{2,Int}[]
    for k = 1:n
        push!(ret, (n-k+1,k))
    end
    ret
end