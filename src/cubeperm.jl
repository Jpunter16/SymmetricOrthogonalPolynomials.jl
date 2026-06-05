function build_A(n)
    A = zeros(Int64, 3*n, 3*n)
    A[1, 2*n] = 1
    for i = 2:2*n-1
        if i % 2 == 0
            A[i, 2*n - i] = 1
        else
            A[i, 2*n + div(i, 2)] = 1
        end
    end
    for i=2*n:3*n-1
        A[i,1+2(i-2*n)]=1
    end
    A[end,end]=1
    A=sparse(A)
    dropzeros(A)
    return A

end

#build the representation applied to the permutation (1) (2 3):

function build_B(n)
    B = zeros(Int64, 3*n, 3*n)
    B[1,1]=1
    for i=2:2n-1
        if i % 2 == 0
            B[i, i+1] = 1
        else
            B[i, i-1] = 1
        end
    end
    for i=2*n:3*n
        B[i,end-(i-2*n)]=1
    end
    B=sparse(B)
    dropzeros(B)
    return B

end

#get the matrix that orthogonalises the representation by blocks
function get_Q(n)
    rep=Representation{SparseMatrixCSC{Float64, Int64}}(SparseMatrixCSC{Float64, Int64}[build_A(n),build_B(n)])
    lam,Q=blockdiagonalize(rep)
    return(Q)
end