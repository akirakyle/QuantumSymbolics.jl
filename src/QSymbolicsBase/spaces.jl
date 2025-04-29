abstract type AbstractHilbertSpace end
abstract type AbstractKetSpace <: AbstractHilbertSpace end
abstract type AbstractBraSpace <: AbstractHilbertSpace end
abstract type AbstractOperatorSpace end
abstract type AbstractSuperOperatorSpace end
const StateSpace = Union{KetSpace, BraSpace}
const LinearSpace = Union{OperatorSpace, SuperOperatorSpace}

struct KetSpace <: AbstractKetSpace
    shape::Vector{Int}
end
struct BraSpace <: AbstractBraSpace
    shape::Vector{Int}
end
Base.:(==)(s1::StateSpace, s2::StateSpace) = a.shape == b.shape
shape(s::StateSpace) = s.shape
dimension(s::StateSpace) = prod(s.shape)

struct OperatorSpace{T<:AbstractHilbertSpace} <: AbstractOperatorSpace
    left::T
    right::T
end
struct SuperOperatorSpace{T<:AbstractHilbertSpace} <: AbstractSuperOperatorSpace
    space_l::T
    space_r::T
end
Base.:(==)(s1::LinearSpace, s2::LinearSpace) = a.space_l == b.space_l && a.space_r == b.space_r
space_l(s::LinearSpace) = s.left
space_r(s::LinearSpace) = s.right
dimension(s::OperatorSpace) = dimension(s.left)*dimension(s.right)
dimension(s::SuperOperatorSpace) = (dimension(s.left)*dimension(s.right))^2

# not correct for oscillators or rotors, need to figure out better scheme

tensor(a::T, b::T) where T<:StateSpace = T([shape(a); shape(b)])
tensor(a::T, b::T) where T<:LinearSpace = T(tensor(space_l(a), space_l(b)), tensor(space_r(a), space_r(b)))

dagger(s::KetSpace) = BraSpace(s.shape)
dagger(s::BraSpace) = KetSpace(s.shape)
dagger(s::T) where T<:LinearSpace = T(dagger(space_l(s)), dagger(space_r(s)))

TrivialSpace() = KetSpace([1])
QubitSpace() = KetSpace([2])
#QubitSpace{T}() where T<:AbstractHilbertSpace = T([2])
#NLevelSpace{T}(N) where T<:AbstractHilbertSpace = T([N])
#OscillatorSpace{T}() where T<:AbstractHilbertSpace = T([0])
#RotorSpace{T}() where T<:AbstractHilbertSpace = T([-1])

TrivialSpace{OperatorSpace{T}}() where T = OperatorSpace(TrivialSpace{T}(), TrivialSpace{T}())
QubitSpace{OperatorSpace{T}}() where T = OperatorSpace(QubitSpace{T}(), QubitSpace{T}())


"""
Exception that should be raised for an illegal algebraic operation.
"""
mutable struct IncompatibleSpaces <: Exception end

const SPACES_CHECK = Ref(true)

"""
    @compatiblespaces

Macro to skip checks for compatible spaces. Useful for `*`, `expect` and similar
functions.
"""
macro compatiblespaces(ex)
    return quote
        SPACES_CHECK[] = false
        local val = $(esc(ex))
        SPACES_CHECK[] = true
        val
    end
end

"""
    check_addible(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects are not addible as
determined by `addible(a, b)`.  Disabled by use of [`@compatiblebases`](@ref)
anywhere further up in the call stack.
"""
function check_addible(a, b)
    if SPACES_CHECK[] && !addible(a, b)
        throw(IncompatibleSpaces())
    end
    return _add_space(a,b)
end

"""
    check_multiplicable(a, b)

Throw an [`IncompatibleBases`](@ref) error if the objects are not multiplicable
as determined by `multiplicable(a, b)`.  Disabled by use of
[`@compatiblebases`](@ref) anywhere further up in the call stack.
"""
function check_multiplicable(a, b)
    if SPACES_CHECK[] && !multiplicable(a, b)
        throw(IncompatibleSpaces())
    end
    return _mul_space(a,b)
end

const AbstractSpace = Union{AbstractHilbertSpace, AbstractOperatorSpace, AbstractSuperOperatorSpace}

adible(a::AbstractSpace, b::AbstractSpace) = (a == b)
addible(a::BasicQSymbolic, b::BasicQSymbolic) = addible(space(a), basis(b))
_add_space(a::AbstractSpace, b::AbstractSpace) = a
_add_space(a::BasicQSymbolic, b::BasicQSymbolic) = space(a)

multiplicible(a::AbstractSpace, b::AbstractSpace) = false
multiplicible(a::BraSpace, b::KetSpace) = a == b
multiplicible(a::KetSpace, b::BraSpace) = true
multiplicible(a::OperatorSpace, b::KetSpace) = space_r(a) == b
multiplicible(a::BraSpace, b::OperatorSpace) = a == space_l(b)
multiplicible(a::OperatorSpace, b::OperatorSpace) = space_r(a) == space_l(b)
multiplicible(a::SuperOperatorSpace, b::OperatorSpace) = space_l(b) == space_r(b) && space_r(a) == space_l(b)
multiplicible(a::OperatorSpace, b::SuperOperatorSpace) = space_l(a) == space_r(a) && space_r(a) == space_l(b)
multiplicible(a::SuperOperatorSpace, b::SuperOperatorSpace) = space_r(a) == space_l(b)
multiplicible(a::BasicQSymbolic, b::BasicQSymbolic) = multiplicible(space(a), space(b))

_mul_space(a::BraSpace, b::KetSpace) = TrivialSpace()
_mul_space(a::KetSpace, b::BraSpace) = OperatorSpace(a,dagger(b))
_mul_space(a::OperatorSpace, b::KetSpace) = b
_mul_space(a::BraSpace, b::OperatorSpace) = a
_mul_space(a::OperatorSpace, b::OperatorSpace) = OperatorSpace(space_l(a), space_r(b))
_mul_space(a::SuperOperatorSpace, b::OperatorSpace) = b
_mul_space(a::OperatorSpace, b::SuperOperatorSpace) = a
_mul_space(a::SuperOperatorSpace, b::SuperOperatorSpace) = SuperOperatorSpace(space_l(a), space_r(b))
_mul_space(a::BasicQSymbolic, b::BasicQSymbolic) = mul_space(space(a), space(b))
