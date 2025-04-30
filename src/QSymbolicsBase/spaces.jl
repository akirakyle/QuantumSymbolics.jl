abstract type AbstractHilbertSpace end
abstract type AbstractOperatorSpace end

struct TrivialSpace <: AbstractHilbertSpace end
struct NLevelSpace{T<: Union{Int, Symbol}} <: AbstractHilbertSpace
    N::T
end
Base.:(==)(a::NLevelSpace, b::NLevelSpace) = a.N == b.N
struct OscillatorSpace <: AbstractHilbertSpace end
struct RotorSpace <: AbstractHilbertSpace end

struct TrivialOpSpace <: AbstractOperatorSpace end
struct NLevelOpSpace{T<: Union{Int, Symbol}} <: AbstractOperatorSpace
    N::T
end
Base.:(==)(a::NLevelOpSpace, b::NLevelOpSpace) = a.N == b.N
struct OscillatorOpSpace <: AbstractOperatorSpace end
struct RotorOpSpace <: AbstractOperatorSpace end

struct LinearSpace{S<:Union{AbstractHilbertSpace, AbstractOperatorSpace}}
    left::Vector{S}
    right::Vector{S}
end
#LinearSpace{T}(hs::S) where {T<:AbstractKet,S} = LinearSpace{T,S}(hs, TrivialSpace())
Base.:(==)(a::LinearSpace, b::LinearSpace) = a.left == b.left && a.right == b.right
#isketspace(s::OperatorSpace) = prod(s.right) == 1
#isbraspace(s::OperatorSpace) = prod(s.left) == 1
#isoperatorspace(s::OperatorSpace) = !(isketspace(s) || isbraspace(s))
#space_l(s::LinearSpace) = s.left
#space_r(s::LinearSpace) = s.right
#dimension(s::OperatorSpace) = dimension(s.left)*dimension(s.right)
#dimension(s::SuperOperatorSpace) = (dimension(s.left)*dimension(s.right))^2
#shape(s::StateSpace) = s.shape
#dimension(s::StateSpace) = prod(s.shape)

tensor(a::LinearSpace{<:AbstractHilbertSpace}, b::LinearSpace{<:AbstractHilbertSpace}) where T = LinearSpace{T}([a.left; b.left], [a.right; b.right])
tensor(a::LinearSpace{<:AbstractOperatorSpace}, b::LinearSpace{<:AbstractOperatorSpace}) where T = LinearSpace{T}([a.left; b.left], [a.right; b.right])
dagger(s::LinearSpace) = LinearSpace(s.right, s.left)

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

adible(a::OperatorSpace, b::OperatorSpace) = (a == b)
addible(a::BasicQSymbolic, b::BasicQSymbolic) = addible(space(a), basis(b))
_add_space(a::OperatorSpace, b::OperatorSpace) = a
_add_space(a::BasicQSymbolic, b::BasicQSymbolic) = space(a)

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

multiplicible(a::LinearSpace, b::LinearSpace) = a.right == b.left
multiplicible(a::BasicQSymbolic, b::BasicQSymbolic) = multiplicible(space(a), space(b))
_mul_space(a::LinearSpace, b::LinearSpace) = LinearSpace(a.left,b.right)
_mul_space(a::BasicQSymbolic, b::BasicQSymbolic) = _mul_space(space(a), space(b))

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
