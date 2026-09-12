import BlockArrays: blocksize

struct S2Invariant{T,B} <: MultivariateOrthogonalPolynomial{2,T}
    basis::B
end
S2Invariant(B::AbstractQuasiMatrix{T}) where {T} = S2Invariant{T,typeof(B)}(B)

S2Invariant() = S2Invariant(Legendre())

function productCombination(a::Vector{Int}, b::Vector{Int})
    if length(a)!=2
        error("incorrrect length a")
    elseif length(b)!=2
        error("incorrrect length b")
    end
    v=[[[a[1], b[1]], [a[2], b[2]]], [[a[1], b[2]], [a[2], b[1]]]]
    map(t->map(p->p .+ 1, t), v)
end

function productCombination(a::Partition2, b::Partition2)
    productCombination(a.p, b.p)
end

function massPlusStiff(p, q, P)
    p=map(t->t+1, p)
    q=map(t->t+1, q)
    ###ONLY FOR (1-X^2)C^3/2 basis or C^-1/2
    stiff_mat=P'*diff(P, 2)
    mass_mat=(P'*P)
    stiff_mat[p[1], q[1]]*mass_mat[p[2], q[2]]+mass_mat[p[1], q[1]]*stiff_mat[p[2], q[2]]
end


function getLaplacianS2InvariantBasis(Q::S2Invariant, n::Int)
    a = blockedrange(floor.(Int, ((0:n) ./ 2 .+ 1)))
    P=Q.basis
    Δ=BlockedMatrix(Zeros((a, a)))
    Δ_dim=blocksize(Δ)

    for i = 1:Δ_dim[1]
        part_row_colection=Partition_2_parts(i-1)
        for j = 1:Δ_dim[2]
            part_col_colection=Partition_2_parts(j-1)
            for ib in eachindex(part_row_colection)
                part_row=part_row_colection[ib]
                for jb in eachindex(part_col_colection)
                    part_col=part_col_colection[jb]
                    #prod_comb=productCombination(part_row,part_col)
                    #prod_comb=[collect(zip(part_row.p, p)) for p in permutations(part_col)]
                    norm_row = sqrt(2 + 2*(part_row[1]==part_row[2]))
                    norm_col = sqrt(2 + 2*(part_col[1]==part_col[2]))
                    view(Δ, Block(i, j))[ib, jb] =
                        2/(norm_row*norm_col) * (
                            massPlusStiff(part_row, part_col, P) +
                            massPlusStiff(part_row, reverse(part_col), P) +
                            massPlusStiff(reverse(part_row), part_col, P) +
                            massPlusStiff(reverse(part_row), reverse(part_col), P)
                        )
                end
            end
        end

    end
    Δ
end
