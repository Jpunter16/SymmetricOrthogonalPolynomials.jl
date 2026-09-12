###
# ℤ_2 symmetries
###


function reflection_trivialfilter(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 1:2:N, k = 1:2:K
        ret[Block(K)[k]] .= true
    end
    ret
end

function reflection_signfilter(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 1:2:N, k = 2:2:K
        ret[Block(K)[k]] .= true
    end
    ret
end


function reflection_tsfilter(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 2:2:N, k = 2:2:K
        ret[Block(K)[k]] .= true
    end
    ret
end

function reflection_stfilter(N)
    ret = BlockedArray{Bool}(undef, 1:N)
    ret .= false
    for K = 2:2:N, k = 1:2:K
        ret[Block(K)[k]] .= true
    end
    ret
end
