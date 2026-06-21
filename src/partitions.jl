function Partition_3_parts(n::Integer)
    if n==0
        return [[0,0,0]]
    else
        parts = collect(Combinatorics.partitions(n))
        parts = filter(p -> length(p) <= 3, parts)
        [vcat(p, zeros(Int, 3 - length(p))) for p in parts]
    end
end