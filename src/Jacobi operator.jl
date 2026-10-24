c(n)=n/(2*n+3)

a(n)=0

b(n)=(n+3)/(2*n+3)

function buildPermutation(vals::Vector{Int},mult::Vector{Int})
    reduce(vcat, fill.(vals, mult))
end

function JacobiOperatorE1(α::Vector{Int})
    vals = sort(unique(α); rev = true)
    mult = [count(==(v), α) for v in vals]

    buckets = Dict{Vector{Int},Float64}()

    for i in eachindex(vals)
        if vals[i]>=1
            new_ind=vals[i]-1
            new_mult=copy(mult)
            new_mult[i]-=1
            new_part=buildPermutation(vals,new_mult)
            append!(new_part,new_ind)
            sort!(new_part; rev=true)
            buckets[new_part]=get(buckets, new_part, 0 ) + c(vals[i])*mult[i]
        end
        new_ind=vals[i]+1
        new_mult=copy(mult)
        new_mult[i]-=1
        new_part=buildPermutation(vals,new_mult)
        append!(new_part,new_ind)
        sort!(new_part; rev=true)
        buckets[new_part]=get(buckets, new_part, 0 ) + b(vals[i])*mult[i]
    end
    buckets
end
