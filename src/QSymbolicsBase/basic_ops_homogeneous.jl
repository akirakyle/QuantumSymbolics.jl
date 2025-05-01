function _set_add_dict_kv!(d, k, v)
    ak = get(d, k, nothing)
    if ak !== nothing
        v = ak + v
    end
    if iszero(v)
        delete!(d, k)
    else
        d[k] = v
    end
end

function _set_add_dict!(d, x)
    if isadd(x)
        for (k, v) in x.dict
            _set_add_dict_kv!(d,k,v)
        end
    else
        k,v = get_coeff_and_op(x)
        _set_add_dict_kv!(d,k,v)
    end
end

function +(a::BasicQSymbolic{T}, b::BasicQSymbolic{S}) where {T,S}
    check_addible(a,b)

    if isadd(a)
        d = Dict(a.dict)
        _set_add_dict!(d, b)
    elseif isadd(b)
        d = Dict(b.dict)
        _set_add_dict!(d, a)
    else
        d = Dict{BasicQSymbolic, Symbolic}()
        _set_add_dict!(d, a)
        _set_add_dict!(d, b)
    end
    length(d) == 0 && return zero(a)
    add = Add{T}(dict=d, basis_l=basis_l(a), basis_r=basis_r(a))
    metadata(add)[:hermitian] = ishermitian(a) && ishermitian(a)
end

Base.:(-)(x::BasicQSymbolic) = (-1)*x
Base.:(-)(x::BasicQSymbolic,y::BasicQSymbolic) = x + (-y)

function Base.:(*)(c::T, x::BasicQSymbolic{S}) where {T<:Union{Number, Symbolic{<:Number}}, S}
    if (isa(c, Number) && iszero(c))
        zero(x)
    elseif (isa(c, Number) && isone(c))
        x
    else
        bases = (:basis_l => basis_l(x), :basis_r => basis_r(x))
        op = @match x begin
            (Sym(_) || Add(_)) => Mul{S}(coeff=c, terms=[x], bases...)
            Mul(_) => Mul{S}(coeff=c*x.coeff, terms=x.terms, bases...)
            Tensor(_) => Tensor{S}(coeff=c*x.coeff, terms=x.terms, bases...)
            Sum(_) => Sum{S}(coeff=c*x.coeff, terms=x.terms, bases...)
        end
        metadata(op)[:hermitian] = _isreal(c) && ishermitian(x)
        op
    end
end

Base.:(*)(x::BasicQSymbolic, c::T) where {T<:Union{Number, Symbolic{<:Number}}} = c*x
Base.:(/)(x::BasicQSymbolic, c::T) where {T<:Union{Number, Symbolic{<:Number}}} = iszero(c) ? throw(DomainError(c,"cannot divide QSymbolics expressions by zero")) : (1/c)*x

_mul_type(A::Type{<:AbstractBra}, B::Type{<:AbstractKet}) = Number
_mul_type(A::Type{<:AbstractKet}, B::Type{<:AbstractBra}) = AbstractOperator
_mul_type(A::Type{<:AbstractKet}, B::Type{<:AbstractBra}) = AbstractOperator

function Base.:(*)(a::BasicQSymbolic{T}, b::BasicQSymbolic{S})
    U = _mul_type(T,S)
    check_multiplicable(a,b)
    if iszero(a) || iszero(b)
        return SZero(basis_l(a), basis_r(b))
    end
    bases = (:basis_l => basis_l(a), :basis_r => basis_r(b))
    if ismul(a) && ismul(b)
        op = Mul{U}(coeff=a.coeff*b.coeff, terms=[a.terms; b.terms], bases...)
    elseif ismul(a)
        k,v = get_coeff_and_op(b)
        op = Mul{U}(coeff=k*a.coeff, terms=[a.terms; v], bases...)
    elseif ismul(b)
        k,v = get_coeff_and_op(a)
        op = Mul{U}(coeff=k*b.coeff, terms=[v; b.terms], bases...)
    elseif (istensor(a) && istensor(b)) || (issum(a) && issum(b))
        op = Tensor{U}(coeff=a.coeff*b.coeff, terms=[x*y for (x,y) in zip(a.terms, b.terms)], bases...)
    else
        ka,va = get_coeff_and_op(a)
        kb,vb = get_coeff_and_op(b)
        op = Mul{U}(coeff=ka*kb, terms=[va,vb], bases...)
    end
    metadata(op)[:unitary] = isunitary(a) && isunitary(b)
    op
end

function ⊗(a::BasicQSymbolic, b::BasicQSymbolic)
    s = tensor(space(a), space(b))
    if iszero(a) || iszero(b)
        return SZero(s)
    end
    args = (:hermitian => ishermitian(a) && ishermitian(b),
            :unitary => isunitary(a) && isunitary(b), :space => s)
    if istensor(a) && istensor(b)
        Tensor(coeff=a.coeff*b.coeff, terms=[a.terms; b.terms], args...)
    elseif istensor(a)
        k,v = get_coeff_and_op(b)
        Tensor(coeff=k*a.coeff*b, terms=[a.terms; v], args...)
    elseif istensor(b)
        k,v = get_coeff_and_op(a)
        Mul(coeff=k*b.coeff, terms=[v; b.terms], args...)
    elseif (istensor(a) && istensor(b)) || (issum(a) && issum(b))
        Tensor(coeff=a.coeff*b.coeff, terms=[x*y for (x,y) in zip(a.terms, b.terms)], args...)
    else
        ka,va = get_coeff_and_op(a)
        kb,vb = get_coeff_and_op(b)
        Tensor(coeff=ka*kb, terms=[va,vb], args...)
    end
end

function kraus(xs::SymQO...)
    xs = collect(xs)
    f
end

Base.:(*)(sop::KrausRepr, op::Symbolic{AbstractOperator}) = (+)((i*op*dagger(i) for i in sop.krausops)...)
Base.:(*)(sop::KrausRepr, k::Symbolic{AbstractKet}) = (+)((i*SProjector(k)*dagger(i) for i in sop.krausops)...)
Base.:(*)(sop::KrausRepr, k::SZeroOperator) = SZeroOperator()
Base.:(*)(sop::KrausRepr, k::SZeroKet) = SZeroOperator()

Base.:(*)(sop::Symbolic{AbstractSuperOperator}, op::Symbolic{AbstractOperator}) = SSuperOpApply(sop,op)
Base.:(*)(sop::Symbolic{AbstractSuperOperator}, op::SZeroOperator) = SZeroOperator()
Base.:(*)(sop::Symbolic{AbstractSuperOperator}, k::Symbolic{AbstractKet}) = SSuperOpApply(sop,SProjector(k))
Base.:(*)(sop::Symbolic{AbstractSuperOperator}, k::SZeroKet) = SZeroOperator()

basis(x::SVec) = (⊗)(fill(basis(x.op), length(basis(x.op)))...)
Base.show(io::IO, x::SVec) = print(io, "|$(x.op)⟩⟩")

vec(x::Symbolic{AbstractOperator}) = SVec(x)
vec(x::SScaled{AbstractOperator}) = x.coeff*vec(x.obj)
vec(x::SAdd{AbstractOperator}) = (+)((vec(i) for i in arguments(x))...)
