using NumericalRepresentationTheory

abstract type PartitionGeneralized end

struct Partition3 <:PartitionGeneralized
    p::Vector{Int}
    function Partition3(p)
        if length(p)!=3
            error("vector should be of length 3")
        end
        if !issorted(p; lt=Base.:>)
            error("input vector $p should be sorted")
        end
        if !all(x -> x >= 0, p)
            error("input vector $p should be nonnegative")
        end
        new(p)
    end
end


function Partition_3_parts(n::Int)
    if n==0
        return [Partition3([0,0,0])]
    else
        parts = collect(Combinatorics.partitions(n))
        parts = filter(p -> length(p) <= 3, parts)
        [Partition3(vcat(p, zeros(Int, 3 - length(p)))) for p in parts]
    end
end

function getindex(p::Partition3,n::Int)
    p.p[n]
end

struct Partition2 <:PartitionGeneralized
    p::Vector{Int}
    function Partition2(p)
        if length(p)!=2
            error("vector should be of length 3")
        end
        if !issorted(p; lt=Base.:>)
            error("input vector $p should be sorted")
        end
        if !all(x -> x >= 0, p)
            error("input vector $p should be all nonnegative")
        end
        new(p)
    end
end

function Partition_2_parts(n::Int)
    if n==0
        return [Partition2([0,0])]
    else
        parts = collect(Combinatorics.partitions(n))
        parts = filter(p -> length(p) <= 2, parts)
        [Partition2(vcat(p, zeros(Int, 2 - length(p)))) for p in parts]
    end
end

function getindex(p::Partition2,n::Int)
    p.p[n]
end

function Base.convert(::Type{NumericalRepresentationTheory.Partition}, y::Partition2)
    NumericalRepresentationTheory.Partition(filter(t->t !=0, y.p))
end
function Base.convert(::Type{NumericalRepresentationTheory.Partition}, y::Partition3)
    NumericalRepresentationTheory.Partition(filter(t->t !=0, y.p))
end