abstract type AbstractSpace end
abstract type AbstractHilbertSpace <: AbstractSpace end
abstract type AbstractOperatorSpace <: AbstractSpace end
abstract type AbstractSuperOpSpace <: AbstractSpace end

struct TrivialSpace <: AbstractHilbertSpace end
struct OscillatorSpace <: AbstractHilbertSpace end
struct RotorSpace <: AbstractHilbertSpace end
struct NLevelSpace{T<: Union{Int, Symbol}} <: AbstractHilbertSpace
    N::T
end
Base.:(==)(a::NLevelSpace, b::NLevelSpace) = a.N == b.N
QubitSpace() = NLevelSpace(2)

struct CompositeHilbertSpace{T<:AbstractHilbertSpace} <: AbstractHilbertSpace
    s::Vector{T}
end
Base.:(==)(a::CompositeHilbertSpace, b::CompositeHilbertSpace) = a.s == b.s
space_l(s::QOSpace) = s.left
space_r(s::QOSpace) = s.right

struct QOSpace{T<:Union{AbstactHilbertSpace, AbstractOperatorSpace}} <: AbstractSpace
    left::T
    right::T
end
Base.:(==)(a::QOSpace, b::QOSpace) = a.left == b.left && a.right == b.right
space_l(s::QOSpace) = s.left
space_r(s::QOSpace) = s.right
#dimension(s::OperatorSpace) = dimension(s.left)*dimension(s.right)
#dimension(s::SuperOperatorSpace) = (dimension(s.left)*dimension(s.right))^2
#shape(s::StateSpace) = s.shape
#dimension(s::StateSpace) = prod(s.shape)

tensor(a::QOSpace{<:AbstractHilbertSpace}, b::QOSpace{<:AbstractHilbertSpace}) where T = QOSpace{T}([a.left; b.left], [a.right; b.right])
tensor(a::QOSpace{<:AbstractOperatorSpace}, b::QOSpace{<:AbstractOperatorSpace}) where T = QOSpace{T}([a.left; b.left], [a.right; b.right])
dagger(s::QOSpace) = QOSpace(s.right, s.left)

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

addible(a::BasicQSymbolic, b::BasicQSymbolic) = space(a) == basis(b)
_add_type_space(a::BasicQSymbolic{T}, b::BasicQSymbolic{T}) = (T, space(a))

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
end

multiplicable(a::BasicQSymbolic, b::BasicQSymbolic) = space(a).right == space(b).left
_mul_type_space(a::BasicQSymbolic{<:AbstractKet}, b::BasicQSymbolic{<:AbstractBra}) =
    (AbstractOperator, QOSpace(space(a).left, space(b).right))
_mul_type_space(a::BasicQSymbolic{<:AbstractBra}, b::BasicQSymbolic{<:AbstractKet}) =
    (AbstractOperator, QOSpace(space(a).left, space(b).right))
_mul_type_space(a::BasicQSymbolic{<:AbstractOperator}, b::BasicQSymbolic{<:AbstractKet}) =
    (AbstractKet, QOSpace(space(a).left, space(b).right))

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
end
