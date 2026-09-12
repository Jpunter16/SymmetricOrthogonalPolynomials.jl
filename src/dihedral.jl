
"""
    dihedralQ(n)

gives the orthogonal transformation from tensor product to Dihedral symmetry-respecting basis.
"""

dihedralQ(n) = dihedralQ(Float64, n)
function dihedralQ(T, n)
    Q = zeros(T, n, n)
    if iseven(n)
        for k = 1:2:(n÷2)
            Q[k, k] = 1
            Q[k+1, n-k+1] = 1
        end
        m = last(1:2:(n÷2))+1
        for k = 1:2:(n÷2-1)
            Q[k+m+1, k+1] = 1
            Q[k+m, n-k] = 1
        end
    elseif mod(n, 4) == 1
        for k = 1:(n÷4)
            Q[k, 2k-1] = Q[k, n-2k+2] = 1/sqrt(T(2))
        end
        m = n÷4
        Q[m+1, n÷2+1] = 1
        m = n÷4+1
        for k = 1:(n÷4)
            Q[k+m, 2k-1] = 1/sqrt(T(2))
            Q[k+m, n-2k+2] = -1/sqrt(T(2))
        end
        m = n÷2+1
        for k = 1:(n÷4)
            Q[k+m, 2k] = 1/sqrt(T(2))
            Q[k+m, n-2k+1] = -1/sqrt(T(2))
        end
        m = 3*(n÷4)+1
        for k = 1:(n÷4)
            Q[k+m, 2k] = 1/sqrt(T(2))
            Q[k+m, n-2k+1] = 1/sqrt(T(2))
        end
    elseif mod(n, 4) == 3
        for k = 1:(n÷4+1)
            Q[k, 2k-1] = Q[k, n-2k+2] = 1/sqrt(T(2))
            m = n÷4+1
            Q[k+m, 2k-1] = 1/sqrt(T(2))
            Q[k+m, n-2k+2] = -1/sqrt(T(2))
        end
        m = n÷2+1
        for k = 1:(n÷4)
            Q[k+m, 2k] = 1/sqrt(T(2))
            Q[k+m, n-2k+1] = -1/sqrt(T(2))
        end
        m = 3*(n÷4)+2
        for k = 1:(n÷4)
            Q[k+m, 2k] = 1/sqrt(T(2))
            Q[k+m, n-2k+1] = 1/sqrt(T(2))
        end
        Q[end, n÷2+1] = 1
    end
    Q
end



function dihedral_trivialfilter(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 1:2:N
        ret[Block(K)[1:((K-1)÷4+1)]] .= true
    end
    ret
end

function dihedral_signfilter(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 3:2:N
        ret[Block(K)[(K-(K-2)÷4):K]] .= true
    end
    ret
end

function dihedral_tsfilter(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 3:2:N
        ret[Block(K)[((K-1)÷4+2):(K÷2+1)]] .= true
    end
    ret
end

function dihedral_stfilter(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 5:2:N
        ret[Block(K)[(K÷2+2):(K-(K-2)÷4-1)]] .= true
    end
    ret
end

function dihedral_faithfulfilter1(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 2:2:N, k = 1:2:K
        ret[Block(K)[k]] .= true
    end
    ret
end


function dihedral_faithfulfilter2(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 2:2:N, k = 2:2:K
        ret[Block(K)[k]] .= true
    end
    ret
end
