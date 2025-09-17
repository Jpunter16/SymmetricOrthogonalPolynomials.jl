function cubevectorrepresentation(n)
    m = 3length(lextuples(n))
    S₁ = BlockArray(zeros(Int,m,m), Fill(m÷3,3), Fill(m÷3,3))
    S₁[Block(1,3)] = cubegen1(n)
    S₁[Block(2,2)] = -cubegen1(n)
    S₁[Block(3,1)] = cubegen1(n)
    S₂ = BlockArray(zeros(Int,m,m), Fill(m÷3,3), Fill(m÷3,3))
    S₂[Block(1,1)] = -cubegen2(n)
    S₂[Block(2,3)] = -cubegen2(n)
    S₂[Block(3,2)] = -cubegen2(n)
    S₃ = BlockArray(zeros(Int,m,m), Fill(m÷3,3), Fill(m÷3,3))
    S₃[Block(3,1)] = -cubegen3(n)
    S₃[Block(2,2)] = -cubegen3(n)
    S₃[Block(1,3)] = -cubegen3(n)
    Representation(Matrix{Int}[S₁,S₂,S₃])
end

function cubevector_filter((p, s), N, j...)
    ret = zeros(Bool, 3binomial(N+2, N-1))
    ind = 0
    for n = 1:N
        M = 3sum(1:n)
        if s == iseven(n)
            _cube_filter!(view(ret, ind+1:ind+M), cubevectorrepresentation(n), p, n, j...)
        end
        ind += M
    end
    ret
end