import BlockArrays:blocksize

struct S2Invariant{T,B} <: MultivariateOrthogonalPolynomial{2,T}
    basis::B
end
S2Invariant(B::AbstractQuasiMatrix{T}) where T = S2Invariant{T, typeof(B)}(B)

S2Invariant() = S2Invariant(Legendre())

function productCombination(a::Vector{Int},b::Vector{Int})
    if length(a)!=2
        error("incorrrect length a")
    elseif length(b)!=2
        error("incorrrect length b")
    end
    v=[[[a[1],b[1]],[a[2],b[2]]],[[a[1],b[2]],[a[2],b[1]]]]
    map(t->map(p->p.+1,t),v)
end

function productCombination(a::Partition2,b::Partition2)
    productCombination(a.p,b.p)
end


function getLaplacianS2InvariantBasis(Q::S2Invariant, n::Int)
    a = blockedrange(floor.(Int,((0:n)./2 .+1)) )
    P=Q.basis
    Δ=BlockedMatrix(Zeros((a,a)))
    Δ_dim=blocksize(Δ)

    stiff_mat=(P'*diff(P,2))[1:n+1,1:n+1]
    mass_mat=(P'*P)[1:n+1,1:n+1]
    for i=1:Δ_dim[1]
        part_row_colection=Partition_2_parts(i-1)
        for j=1:Δ_dim[2] 
            part_col_colection=Partition_2_parts(j-1)
            for ib=eachindex(part_row_colection)
                part_row=part_row_colection[ib]
                for  jb=eachindex(part_col_colection) 
                    part_col=part_col_colection[jb]
                    prod_comb=productCombination(part_row,part_col)
                    norm_row = sqrt(2 + 2*(part_row[1]==part_row[2]))
                    norm_col = sqrt(2 + 2*(part_col[1]==part_col[2]))
                    view(Δ,Block(i,j))[ib,jb]=
                    2/(norm_row*norm_col) *
                    (mass_mat[prod_comb[1][1]...]*stiff_mat[prod_comb[1][2]...]
                    + mass_mat[prod_comb[2][1]...]*stiff_mat[prod_comb[2][2]...]
                    + stiff_mat[prod_comb[2][1]...]*mass_mat[prod_comb[2][2]...]
                    + stiff_mat[prod_comb[1][1]...]*mass_mat[prod_comb[1][2]...])
                end
            end
        end
        
    end
    Δ
end