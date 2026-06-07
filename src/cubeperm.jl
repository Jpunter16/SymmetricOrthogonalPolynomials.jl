using SparseArrays, BlockDiagonals

function check_indices(N)
    A = Matrix{Int8}(undef, 0, 3)
    for n=1:N
        iter=1
        for i=reverse(1:n)
            k=0
            for j=reverse(1:(n-i+1)) 
                k+=1
                A=vcat(A, [i-1,j-1,k-1]')
                iter+=1
            end
        end
    end
    collect(eachrow(A))[end-Int(N*(N+1)/2)+1:end,:]
end

check_indices(2)

function make_second_cycle_S3(v)
    vec=copy(v)
    a=vec[2]
    vec[2]=vec[3]
    vec[3]=a
    vec
end


function get_B(N)
    dim=Int(N*(N+1)/2)
    B = Matrix{Int8}(zeros( dim, dim))
    basis_vector_indices=check_indices(N)
    permuted_vector=map(make_second_cycle_S3, basis_vector_indices)
    for i=1:size(B,1)
        pos = findfirst(row -> row == basis_vector_indices[i,:], eachrow(permuted_vector))
        if pos !=Nothing
            B[i, pos]=1
        else
            error("Permutated vector not found")
        end

    end
    B
end

function make_first_cycle_S3(v)
    vec=copy(v)
    a=vec[2]
    vec[2]=vec[1]
    vec[1]=a
    vec
end

function get_A(N)
    dim=Int(N*(N+1)/2)
    A = Matrix{Int8}(zeros( dim, dim))
    basis_vector_indices=check_indices(N)
    permuted_vector=map(make_first_cycle_S3, basis_vector_indices)
    for i=1:size(A,1)
        pos = findfirst(row -> row == basis_vector_indices[i,:], eachrow(permuted_vector))
        if pos !=Nothing
            A[i, pos]=1
        else
            error("Permutated vector not found")
        end

    end
    A
end

function get_Q(n)
    N=1:n
    nrows= Int.(N .*(N .+1) ./2)
    block_diagonals = [Matrix{Float64}(undef, nrows[i], nrows[i]) for i  in 1:n]
    for i=N
        rep= Representation{SparseMatrixCSC{Float64, Int64}}(SparseMatrixCSC{Float64, Int64}[get_A(i),get_B(i)])
        lam,Q=blockdiagonalize(rep)
        block_diagonals[i]=Q
    end
    Block_Q=BlockDiagonal(block_diagonals)
    sparse(Block_Q)
end