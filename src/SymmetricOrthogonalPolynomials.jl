module SymmetricOrthogonalPolynomials
using BlockArrays, NumericalRepresentationTheory, Combinatorics
using MultivariateOrthogonalPolynomials, InfiniteArrays, LazyArrays, DomainSets, StaticArrays, ClassicalOrthogonalPolynomials, BandedMatrices, BlockBandedMatrices, QuasiArrays, LinearAlgebra, ArrayLayouts, ContinuumArrays, SparseArrays
import ContinuumArrays: Basis, grammatrix, @simplify
import BlockArrays: block, blockindex, viewblock
import BlockBandedMatrices: AbstractBandedBlockBandedMatrix, blockbandwidths, subblockbandwidths
using MultivariateOrthogonalPolynomials: MultivariateOrthogonalPolynomial, BlockOneTo
import Base: axes, getindex, size
import FastTransforms:pochhammer
import BandedMatrices: _BandedMatrix

export dihedralQ, dihedral_signfilter, dihedral_trivialfilter, dihedral_tsfilter, dihedral_stfilter, dihedral_faithfulfilter1, dihedral_faithfulfilter2
export reflection_trivialfilter, reflection_signfilter, reflection_tsfilter, reflection_stfilter
export cuberepresentation, cube_filter, lextuples, lextuples2

export Partition_3_parts, getLaplacianS2InvariantBasis, Partition_2_parts

export DihedralInvariant, DihedralWeakLaplacian, S3Invariant, S2Invariant, get_Q

export cubeperm_inds

include("dihedral.jl")
include("reflection.jl")
include("cube.jl")
include("square.jl")

include("cubevector.jl")

include("partitions.jl")

include("dihedralinvariant.jl")
include("S3invariant.jl")
include("S2Invariant.jl")
include("cubeperm.jl")
include("RepTheoryPDEs.jl")



end # module