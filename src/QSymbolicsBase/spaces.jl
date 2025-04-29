
struct OperatorSpace{T} where T <: Union{Int, Symbol}
    left::Vector{T}
    right::Vector{T}
end
Base.:(==)(s1::OperatorSpace, s2::OperatorSpace) = a.left == b.left && a.right == b.right
#space_l(s::LinearSpace) = s.left
#space_r(s::LinearSpace) = s.right
#dimension(s::OperatorSpace) = dimension(s.left)*dimension(s.right)
#dimension(s::SuperOperatorSpace) = (dimension(s.left)*dimension(s.right))^2
#shape(s::StateSpace) = s.shape
#dimension(s::StateSpace) = prod(s.shape)

isketspace(s::OperatorSpace) = prod(s.right) == 1
isbraspace(s::OperatorSpace) = prod(s.left) == 1
isoperatorspace(s::OperatorSpace) = !(isketspace(s) || isbraspace(s))

# not correct for oscillators or rotors, need to figure out better scheme
# use symbol for these

tensor(a::T, b::T) where T<:OperatorSpace = T([a.left; b.left], [a.right; b.right])
dagger(s::OperatorSpace) = OperatorSpace(s.right, s.left)

function _make_space{T}(N) where T
    if T <: AbstractKet
        OperatorSpace([N], [1])
    elseif T <: AbstractBra
        OperatorSpace([1], [N])
    elseif T <: AbstractOperator
        OperatorSpace([N], [N])
    else
        throw(ArgumentError())
    end
end

TrivialSpace() = OperatorSpace([1], [1])
QubitSpace{T}() where T = _make_space{T}(2)
NLevelSpace{T}(N::Int) = _make_space{T}(N)
OscillatorSpace{T}() where T = _make_space{T}(:oscillator)
RotorSpace{T}() where T = _make_space{T}(:rotor)

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

multiplicible(a::OperatorSpace, b::OperatorSpace) = a.right == b.left
multiplicible(a::BasicQSymbolic, b::BasicQSymbolic) = multiplicible(space(a), space(b))
_mul_space(a::OperatorSpace, b::OperatorSpace) = OperatorSpace(a.left,b.right)
_mul_space(a::BasicQSymbolic, b::BasicQSymbolic) = mul_space(space(a), space(b))

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
