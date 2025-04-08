
##
# TensorBasis, SumBasis
##

"""
    CompositeBasis(b1, b2...)

Basis for composite Hilbert spaces.

Stores the subbases in a vector and creates the shape vector directly from the
dimensions of these subbases. Instead of creating a CompositeBasis directly,
`tensor(b1, b2...)` or `b1 ⊗ b2 ⊗ …` should be used.
"""
struct CompositeBasis{S<:Integer,B<:Basis} <: Basis
    shape::Vector{S}
    bases::Vector{B}
end
CompositeBasis(bases) = CompositeBasis([length(b) for b in bases], bases)
CompositeBasis(bases::Basis...) = CompositeBasis([bases...])
CompositeBasis(bases::Tuple) = CompositeBasis([bases...])

Base.:(==)(b1::CompositeBasis, b2::CompositeBasis) = all(((i, j),) -> i == j, zip(b1.bases, b2.bases))
Base.length(b::CompositeBasis) = prod(b.shape)
Base.size(b::CompositeBasis) = b.shape
Base.getindex(b::CompositeBasis, i) = getindex(b.bases, i)

"""
    tensor(x::Basis, y::Basis, z::Basis...)

Create a [`CompositeBasis`](@ref) from the given bases.

Any given CompositeBasis is expanded so that the resulting CompositeBasis never
contains another CompositeBasis.
"""
tensor(b1::Basis, b2::Basis) = CompositeBasis([length(b1), length(b2)], [b1, b2])
tensor(b1::CompositeBasis, b2::CompositeBasis) = CompositeBasis([b1.shape; b2.shape], [b1.bases; b2.bases])
tensor(b1::CompositeBasis, b2::Basis) = CompositeBasis([b1.shape; length(b2)], [b1.bases; b2])
tensor(b1::Basis, b2::CompositeBasis) = CompositeBasis([length(b1); b2.shape], [b1; b2.bases])
tensor(bases::Basis...) = reduce(tensor, bases)
tensor(basis::Basis) = basis

function Base.:^(b::Basis, N::Integer)
    if N < 1
        throw(ArgumentError("Power of a basis is only defined for positive integers."))
    end
    tensor([b for i=1:N]...)
end

"""
    SumBasis(b1, b2...)

Similar to [`CompositeBasis`](@ref) but for the [`directsum`](@ref) (⊕)
"""
struct SumBasis{S<:Integer,B<:Basis} <: Basis
    shape::Vector{S}
    bases::Vector{B}
end
SumBasis(bases) = SumBasis([length(b) for b in bases], bases)
SumBasis(bases::Basis...) = SumBasis([bases...])
SumBasis(bases::Tuple) = SumBasis([bases...])

Base.:(==)(b1::SumBasis, b2::SumBasis) = all(((i, j),) -> i == j, zip(b1.bases, b2.bases))
Base.length(b::SumBasis) = sum(b.shape)
Base.getindex(b::SumBasis, i) = getindex(b.bases, i)

"""
    directsum(b1::Basis, b2::Basis)

Construct the [`SumBasis`](@ref) out of two sub-bases.
"""
directsum(b1::Basis, b2::Basis) = SumBasis([length(b1), length(b2)], [b1, b2])
directsum(b1::SumBasis, b2::SumBasis) = SumBasis([b1.shape, b2.shape], [b1.bases; b2.bases])
directsum(b1::SumBasis, b2::Basis) = SumBasis([b1.shape; length(b2)], [b1.bases; b2])
directsum(b1::Basis, b2::SumBasis) = SumBasis([length(b1); b2.shape], [b1; b2.bases])
directsum(bases::Basis...) = reduce(directsum, bases)
directsum(basis::Basis) = basis

embed(b::SumBasis, indices, ops) = embed(b, b, indices, ops)

##
# Common bases
##

"""
    FockBasis(N,offset=0)

Basis for a Fock space where `N` specifies a cutoff, i.e. what the highest
included fock state is. Similarly, the `offset` defines the lowest included fock
state (default is 0). Note that the dimension of this basis is `N+1-offset`.
The [`cutoff`](@ref) and [`offset`](@ref) functions can be used to obtain the
respective properties of a given `FockBasis`.
"""
struct FockBasis{T<:Integer} <: Basis
    N::T
    offset::T
    function FockBasis(N::T,offset::T=0) where T
        if N < 0 || offset < 0 || N <= offset
            throw(DimensionMismatch())
        end
        new{T}(N, offset)
    end
end

Base.:(==)(b1::FockBasis, b2::FockBasis) = (b1.N==b2.N && b1.offset==b2.offset)
Base.length(b::FockBasis) = b.N - b.offset + 1

"""
    cutoff(b::FockBasis)

Return the fock cutoff of the given fock basis.

See [`FockBasis`](@ref).
"""
cutoff(b::FockBasis) = b.N

"""
    offset(b::FockBasis)

Return the offset of the given fock basis.

See [`FockBasis`](@ref).
"""
offset(b::FockBasis) = b.offset


"""
    NLevelBasis(N)

Basis for a system consisting of N states.
"""
struct NLevelBasis{T<:Integer} <: Basis
    N::T
    function NLevelBasis(N::T) where T
        if N < 1
            throw(DimensionMismatch())
        end
        new{T}(N)
    end
end

Base.:(==)(b1::NLevelBasis, b2::NLevelBasis) = b1.N == b2.N
Base.length(b::NLevelBasis) = b.N

"""
    SpinBasis(n)

Basis for spin-n particles.

The basis can be created for arbitrary spin numbers by using a rational number,
e.g. `SpinBasis(3//2)`. The Pauli operators are defined for all possible spin
numbers. The [`spinnumber`](@ref) function can be used to get the spin number
for a `SpinBasis`.
"""
struct SpinBasis{T<:Integer} <: Basis
    spinnumber::Rational{T}
    function SpinBasis(spinnumber::Rational{T}) where T
        n = numerator(spinnumber)
        d = denominator(spinnumber)
        d==2 || d==1 || throw(ArgumentError("Can only construct integer or half-integer spin basis"))
        n >= 0 || throw(ArgumentError("Can only construct positive spin basis"))
        N = numerator(spinnumber*2 + 1)
        new{T}(spinnumber)
    end
end
SpinBasis(spinnumber) = SpinBasis(convert(Rational{Int}, spinnumber))

Base.:(==)(b1::SpinBasis, b2::SpinBasis) = b1.spinnumber==b2.spinnumber
Base.length(b::SpinBasis) = numerator(b.spinnumber*2 + 1)

"""
    spinnumber(b::SpinBasis)

Return the spin number of the given spin basis.

See [`SpinBasis`](@ref).
"""
spinnumber(b::SpinBasis) = b.spinnumber


##
# Operator Bases
##

"""
    KetBraBasis(BL,BR)

The "Ket-Bra" operator basis is the standard representation for the left and
right bases of superoperators. This basis is formed by "vec'ing" the
outer-product "Ket-Bra" basis for an operator with a left Bra basis and right
Ket basis which practically means flipping the Bra to a Ket. The operator itself
is then represented as a "Super-Bra" in this basis and corresponds to
column-stacking its matrix.
"""
struct KetBraBasis{BL<:Basis, BR<:Basis} <: Basis
    left::BL
    right::BR
end
KetBraBasis(b::Basis) = KetBraBasis(b,b)
basis_l(b::KetBraBasis) = b.left
basis_r(b::KetBraBasis) = b.right
Base.:(==)(b1::KetBraBasis, b2::KetBraBasis) = (b1.left == b2.left && b1.right == b2.right)
Base.length(b::KetBraBasis) = length(b.left)*length(b.right)
Base.size(b::KetBraBasis) = (length(b.left), length(b.right))

struct ChoiRefSysBasis{B<:Basis} <: Basis
    basis::B
end
Base.:(==)(b1::ChoiRefSysBasis, b2::ChoiRefSysBasis) = (b1.basis == b2.basis)
Base.length(b::ChoiRefSysBasis) = length(b.basis)
Base.size(b::ChoiRefSysBasis) = (length(b.basis),)

struct ChoiOutSysBasis{B<:Basis} <: Basis
    basis::B
end
Base.:(==)(b1::ChoiOutSysBasis, b2::ChoiOutSysBasis) = (b1.basis == b2.basis)
Base.length(b::ChoiOutSysBasis) = length(b.basis)
Base.size(b::ChoiOutSysBasis) = (length(b.basis),)


"""
    _PauliBasis()

Pauli operator basis consisting of the Pauli matrices I, Z, X, Y, in that order.
"""
struct _PauliBasis <: Basis end

Base.:(==)(pb1::_PauliBasis, pb2::_PauliBasis) = true
Base.length(b::_PauliBasis) = 4
Base.size(b::_PauliBasis) = (4,)
