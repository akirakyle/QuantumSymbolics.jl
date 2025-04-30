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

function +(a::BasicQSymbolic, b::BasicQSymbolic)
    s = check_addible(a,b)
    args = (:hermitian => ishermitian(a) && ishermitian(a), :unitary => false, :space => s) 
    if isadd(a)
        d = Dict(a.dict)
        _set_add_dict!(d, b)
        return Add(dict=d, args...)
    elseif isadd(b)
        return b + a
    else
        d = Dict{BasicQSymbolic, Symbolic}()
        _set_add_dict!(d, a)
        _set_add_dict!(d, b)
        return Add(dict=d, args...)
    end
end

Base.:(-)(x::BasicQSymbolic) = (-1)*x
Base.:(-)(x::BasicQSymbolic,y::BasicQSymbolic) = x + (-y)

function Base.:(*)(c::T, x::BasicQSymbolic) where {T<:Union{Number, Symbolic{<:Number}}}
    if (isa(c, Number) && iszero(c))
        zero(x)
    elseif (isa(c, Number) && isone(c))
        x
    else
        args = (:hermitian => _isreal(c) && ishermitian(x), :unitary => false, :space => space(x))
        @match x begin
            (Sym(_) || Add(_)) => Mul(coeff=c, terms=[x], args...)
            Mul(_) => Mul(coeff=c*x.coeff, terms=x.terms, args...)
            Tensor(_) => Tensor(coeff=c*x.coeff, terms=x.terms, args...)
            Sum(_) => Sum(coeff=c*x.coeff, terms=x.terms, args...)
        end
    end
end

Base.:(*)(x::BasicQSymbolic, c::T) where {T<:Union{Number, Symbolic{<:Number}}} = c*x
Base.:(/)(x::BasicQSymbolic, c::T) where {T<:Union{Number, Symbolic{<:Number}}} = iszero(c) ? throw(DomainError(c,"cannot divide QSymbolics expressions by zero")) : (1/c)*x

function Base.:(*)(a::BasicQSymbolic, b::BasicQSymbolic)
    s = check_multiplicable(a,b)
    if iszero(a) || iszero(b)
        return SZero(s)
    end
    args = (:hermitian => false, :unitary => isunitary(a) && isunitary(b), :space => s)
    if ismul(a) && ismul(b)
        Mul(coeff=a.coeff*b.coeff, terms=[a.terms; b.terms], args...)
    elseif ismul(a)
        k,v = get_coeff_and_op(b)
        Mul(coeff=k*a.coeff, terms=[a.terms; v], args...)
    elseif ismul(b)
        k,v = get_coeff_and_op(a)
        Mul(coeff=k*b.coeff, terms=[v; b.terms], args...)
    elseif (istensor(a) && istensor(b)) || (issum(a) && issum(b))
        Tensor(coeff=a.coeff*b.coeff, terms=[x*y for (x,y) in zip(a.terms, b.terms)], args...)
    else
        ka,va = get_coeff_and_op(a)
        kb,vb = get_coeff_and_op(b)
        Mul(coeff=ka*kb, terms=[va,vb], args...)
    end
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

kraus(xs::Symbolic{AbstractOperator}...) = KrausRepr(collect(xs))
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
