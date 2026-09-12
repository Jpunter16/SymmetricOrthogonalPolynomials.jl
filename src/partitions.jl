using NumericalRepresentationTheory

abstract type PartitionGeneralized end

function _validate_partition(p, len::Int)
    if length(p)!=len
        error("vector should be of length $len")
    end
    if !issorted(p; lt = Base.:>)
        error("input vector $p should be sorted")
    end
    if !all(x -> x >= 0, p)
        error("input vector $p should be nonnegative")
    end
end

struct Partition3 <: PartitionGeneralized
    p::Vector{Int}
    function Partition3(p)
        _validate_partition(p, 3)
        new(p)
    end
end

struct Partition2 <: PartitionGeneralized
    p::Vector{Int}
    function Partition2(p)
        _validate_partition(p, 2)
        new(p)
    end
end

function Partition_3_parts(n::Int)
    if n==0
        return [Partition3([0, 0, 0]).p]
    else
        parts = collect(Combinatorics.partitions(n))
        parts = filter(p -> length(p) <= 3, parts)
        [Partition3(vcat(p, zeros(Int, 3 - length(p)))).p for p in parts]
    end
end

function Partition_2_parts(n::Int)
    if n==0
        return [Partition2([0, 0]).p]
    else
        parts = collect(Combinatorics.partitions(n))
        parts = filter(p -> length(p) <= 2, parts)
        [Partition2(vcat(p, zeros(Int, 2 - length(p)))).p for p in parts]
    end
end

struct PartitionN{M} <: PartitionGeneralized
    p::Vector{Int}
    function PartitionN{M}(p) where {M}
        _validate_partition(p, M)
        new{M}(p)
    end
end

function Partition_n_parts(n::Int, m::Int)
    if n==0
        return [PartitionN{m}(zeros(Int, m)).p]
    else
        parts = collect(Combinatorics.partitions(n))
        parts = filter(p -> length(p) <= m, parts)
        [PartitionN{m}(vcat(p, zeros(Int, m - length(p)))).p for p in parts]
    end
end

function getindex(p::PartitionGeneralized, n::Int)
    p.p[n]
end

function PartitionVectorToVector(Partitionv::Vector{<:PartitionGeneralized})
    [part.p for part in Partitionv]
end

function Base.convert(
    ::Type{NumericalRepresentationTheory.Partition},
    y::PartitionGeneralized,
)
    NumericalRepresentationTheory.Partition(filter(t->t != 0, y.p))
end

function Base.reverse(p::PartitionGeneralized)
    reverse(p.p)
end

function has_distinct_parts(p::PartitionGeneralized)
    length(p.p) == length(unique(p.p))
end
