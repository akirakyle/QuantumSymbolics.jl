abstract type QuantumObject end
const SymQO = Symbolic{<:QuantumObject}

struct TrivialSpace <: Basis end
struct OscillatorSpace <: Basis end
struct RotorSpace <: Basis end
struct NLevelSpace{T<: Union{Int, Symbol}} <: Basis
    N::T
end

_default_md() = Dict{Symbol, Any}(
    :express_cache => Dict{Tuple{<:AbstractRepresentation,<:AbstractUse},Any}(),
    :hermitian => false,
    :unitary => false)

# possibly move to simplified typing after this issue is resolved
# https://github.com/Roger-luo/Moshi.jl/issues/33 
@data BasicQSymExpr{T} <: Symbolic{T} begin
    struct Term
        metadata::Dict{Symbol,Any} = _default_md()
        space::QOSpace
        f::Function
        arguments::Vector{BasicQSymExpr.Type} = BasicQSymExpr.Type[]
    end
    struct Sym
        metadata::Dict{Symbol,Any} = _default_md()
        space::QOSpace
        name::Symbol
    end
    struct Add
        metadata::Dict{Symbol,Any} = _default_md()
        space::QOSpace
        dict::Dict{BasicQSymExpr.Type, Symbolic}
    end
    struct Mul
        metadata::Dict{Symbol,Any} = _default_md()
        space::QOSpace
        coeff::Symbolic
        terms::Vector{BasicQSymExpr.Type}
    end
    struct Sum
        metadata::Dict{Symbol,Any} = _default_md()
        space::QOSpace
        coeff::Symboli
        terms::Vector{BasicQSymExpr.Type}
    end
    struct Tensor
        metadata::Dict{Symbol,Any} = _default_md()
        space::QOSpace
        coeff::Symbolic
        terms::Vector{BasicQSymExpr.Type}
    end
end

const BasicQSymbolic = BasicQSymExpr.Type
const Term = BasicQSymExpr.Sym
const Sym = BasicQSymExpr.Sym
const Add = BasicQSymExpr.Add
const Mul = BasicQSymExpr.Mul
const Sum = BasicQSymExpr.Sum
const Tensor = BasicQSymExpr.Tensor

@noinline error_on_type() = error("Internal error: unreachable reached!")
@noinline error_sym() = error("QSym doesn't have a operation or arguments!")

@inline function operation(x::BasicQSymbolic)
    @match x begin
        Term(_) => x.f
        Sym(_) => error_sym()
        Add(_) => (+)
        Mul(_) => (*)
        Sum(_) => (⊕)
        Tensor(_) => (⊗)
    end
end
@inline head(x::BasicQSymbolic) = operation(x)

function arguments(x::BasicQSymbolic)
    @match x begin
        Term(_) => x.arguments
        Sym(_) => error_sym()
        Add(_) => [k*v for k,v in x.dict]
        Mul(_) => [x.cterm; x.qterms]
        Sum(_) => [x.coeff*q for q in x.qterms]
        Tensor(_) => [x.coeff*q for q in x.qterms]
    end
end
children(x::BasicQSymbolic) = arguments(x)

function isexpr(x::BasicQSymbolic)
    @match x begin
        Sym(_) => false
        _ => true
    end
end
iscall(x::BasicQSymbolic) = isexpr(x)

@inline isa_SymType(S, x) = x isa BasicQSymbolic ? variant_name(x) == S : false

"""
    issym(x)

Returns `true` if `x` is a `Sym`. If true, `nameof` must be defined
on `x` and must return a `Symbol`.
"""
issym(x) = isa_SymType(:Sym, x)
isadd(x)  = isa_SymType(:Add, x)
ismul(x)  = isa_SymType(:Mul, x)
istensor(x)  = isa_SymType(:Tensor, x)
issum(x)  = isa_SymType(:Sum, x)

metadata(x::BasicQSymbolic) = x.metadata
space(x::BasicQSymbolic) = x.space
ishermitian(x::BasicQSymbolic) = x.metadata[:hermitian]
isunitary(x::BasicQSymbolic) = x.unitary[:unitary]

addible(a::BasicQSymbolic, b::BasicQSymbolic) = addible(warp(a), wrap(b))
multiplicable(a::BasicQSymbolic, b::BasicQSymbolic) = multiplicable(warp(a), wrap(b))

Base.nameof(s::BasicSymbolic) = issym(s) ? s.name : error("This BasicQSymbolic doesn't have a name")
Base.hash(x::BasicQSymbolic, h::UInt) = isexpr(x) ? hash((head(x), arguments(x)), h) :
    hash((typeof(x),x.space,x.name,x.f,x.arguments), h)
maketerm(::BasicQSymbolic, f, a, m) = f(a...)

function Base.isequal(a::BasicQSymbolic{T}, b::BasicQSymbolic{S}) where {T,S}
    a === b && return true
    variant_name(a) === variant_name(b) || return false
    T === S || return false

    @match (a, b) begin
        (Term(_), Term(_)) => a.f === b.f && isequal(a.arguments, b.arguments)
        (Sym(_), Sym(_)) => a.name === b.name
        (Add(_), Add(_)) => isequal(a.dict, b.dict)
        (Mul(_), Mul(_)) || (Sum(_), Sum(_)) || (Tensor(_), Tensor(_)) =>
            isequal(a.coeffs, b.coeffs) && isequal(a.terms, b.terms)
        _ => error_on_type()
    end
end

#Base.isequal(::BasicQSymbolic, ::Symbolic{Complex}) = false
#Base.isequal(::Symbolic{Complex}, ::BasicQSymbolic) = false

Base.iszero(x::BasicQSymbolic) = isqsym(x) && x.f == zero

function get_coeff_and_op(x::BasicQSymbolic)
    @match x begin
        (Term(_) || Sym(_) || Add(_)) => (1, x)
        Mul(_) => (x.coeff, Mul(hermitian=x.hermitian, unitary=x.unitary, space=x.space, terms=x.terms))
        Sum(_) => (x.coeff, Sum(hermitian=x.hermitian, unitary=x.unitary, space=x.space, terms=x.terms))
        Tensor(_) => (x.coeff, Tensor(hermitian=x.hermitian, unitary=x.unitary, space=x.space, terms=x.terms))
    end
end

#"""Wrapped symbolic bra"""
#@symbolic_wrap struct SBra <: AbstractBra
#    x::Symbolic{AbstractBra}
#end
#SBra(name, basis) = SBra(Sym{AbstractBra}(name=name, basis_l=TrivialSpace(), basis_r=basis))
#basis(bra::SBra) = bra.x.basis_r

SBra(name, basis) = Sym{AbstractBra}(name=name, space=QOSpace{AbstractBra}(basis))
SKet(name, basis) = Sym{AbstractKet}(name=name, space=QOSpace{AbstractKet}(basis))
SOperator(name, basis) = Sym{AbstractOperator}(name, space=QOSpace{AbstractOperator}(KetBraBasis(basis)))
SOperator(name, bl, br) = Sym{AbstractOperator}(name, space=QOSpace{AbstractOperator}(KetBraBasis(bl,br)))
SSuperOperator(name, basis) = Sym{AbstractSuperOperator}(name, space=QOSpace{AbstractSuperOperator}(KetBraBasis(KetBraBasis(basis))))

"""
    @bra(name, space)

Define a symbolic bra of type `SBra`.

```jldoctest
julia> @bra b₁ QubitSpace()
⟨b₁|

julia> @bra b₂ OscillatorSpace()
⟨b₂|
```
"""
macro bra(name, basis)
    :($(esc(name)) = SBra($(Expr(:quote, name)), $(basis)))
end

"""
    @ket(name, space)

Define a symbolic ket of type `SKet`.

```jldoctest
julia> @ket k₁ QubitSpace()
|k₁⟩

julia> @ket k₂ OscillatorSpace()
|k₂⟩
```
"""
macro ket(name, basis)
    :($(esc(name)) = SKet($(Expr(:quote, name)), $(basis)))
end

"""
    @op(name, space)

Define a symbolic operator of type `SOperator`.

```jldoctest
julia> @op A QubitSpace()
A

julia> @op B OscillatorSpace()
B
```
"""
macro op(name, basis)
    :($(esc(name)) = SOperator($(Expr(:quote, name)), $(basis)))
end

macro superop(name, basis)
    :($(esc(name)) = SSuperOperator($(Expr(:quote, name)), $(basis)))
end


# can maybe do something like
# Operator(space, :symbol_name)
# gives a wrapped SOperator with linear algebra defined so tha
# everything wrapps to Ket/Bra/Operator typs with the .data field being a Symbolic{{Ket/Bra/Operator}}
