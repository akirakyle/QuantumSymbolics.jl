
Base.show(io::IO, x::SKet) = print(io, "|$(symbollabel(x))⟩")
Base.show(io::IO, x::SBra) = print(io, "⟨$(symbollabel(x))|")
Base.show(io::IO, x::Union{SOperator,SHermitianOperator,SUnitaryOperator,SHermitianUnitaryOperator,SSuperOperator}) = print(io, "$(symbollabel(x))")
Base.show(io::IO, x::SymQObj) = print(io, symbollabel(x)) # fallback that probably is not great

const SScaledKet = SScaled{AbstractKet}
function Base.show(io::IO, x::SScaledKet)
    if x.coeff isa Real
        print(io, "$(x.coeff)$(x.obj)")
    else
        print(io, "($(x.coeff))$(x.obj)")
    end
end
const SScaledOperator = SScaled{AbstractOperator}
function Base.show(io::IO, x::SScaledOperator)
    if x.coeff isa Real
        print(io, "$(x.coeff)$(x.obj)")
    else
        print(io, "($(x.coeff))$(x.obj)")
    end
end
const SScaledBra = SScaled{AbstractBra}
function Base.show(io::IO, x::SScaledBra)
    if x.coeff isa Real
        print(io, "$(x.coeff)$(x.obj)")
    else
        print(io, "($(x.coeff))$(x.obj)")
    end
end
const SAddBra = SAdd{AbstractBra}
function Base.show(io::IO, x::SAddBra)
    ordered_terms = sort([repr(i) for i in arguments(x)])
    print(io, join(ordered_terms,"+")::String) # type assert to help inference
end
const SAddKet = SAdd{AbstractKet}
function Base.show(io::IO, x::SAddKet)
    ordered_terms = sort([repr(i) for i in arguments(x)])
    print(io, join(ordered_terms,"+")::String) # type assert to help inference
end
const SAddOperator = SAdd{AbstractOperator}
function Base.show(io::IO, x::SAddOperator)
    repr_func = x -> x isa STensor ? "("*repr(x)*")" : repr(x)
    ordered_terms = sort([repr_func(i) for i in arguments(x)])
    print(io, join(ordered_terms,"+")::String) # type assert to help inference
end
function Base.show(io::IO, x::SMulOperator)
    str_func = x -> x isa SAdd || x isa STensor ? "("*string(x)*")" : string(x)
    print(io, join(map(str_func, arguments(x)),""))
end
const STensorBra = STensor{AbstractBra}
Base.show(io::IO, x::STensorBra) = print(io, join(map(string, arguments(x)),""))
const STensorKet = STensor{AbstractKet}
Base.show(io::IO, x::STensorKet) = print(io, join(map(string, arguments(x)),""))
const STensorOperator = STensor{AbstractOperator}
function Base.show(io::IO, x::STensorOperator)
    str_func = x -> x isa SAdd ? "("*string(x)*")" : string(x)
    print(io, join(map(str_func, arguments(x)),"⊗"))
end
const STensorSuperOperator = STensor{AbstractSuperOperator}
function Base.show(io::IO, x::STensorSuperOperator)
    str_func = x -> x isa SAdd ? "("*string(x)*")" : string(x)
    print(io, join(map(str_func, arguments(x)),"⊗"))
end
function Base.show(io::IO, x::SApplyKet) 
    str_func = x -> x isa SAdd || x isa STensorOperator ? "("*string(x)*")" : string(x)
    print(io, join(map(str_func, arguments(x)),""))
end
function Base.show(io::IO, x::SApplyBra) 
    str_func = x -> x isa SAdd || x isa STensor ? "("*string(x)*")" : string(x)
    print(io, join(map(str_func, arguments(x)),""))
end
Base.show(io::IO, x::SBraKet) = begin print(io,x.bra); print(io,x.ket) end
Base.show(io::IO, x::SOuterKetBra) = begin print(io, x.ket); print(io, x.bra) end
Base.show(io::IO, x::KrausRepr) = print(io, "𝒦("*join([symbollabel(i) for i in x.krausops], ",")*")")
Base.show(io::IO, x::SSuperOpApply) = print(io, "$(x.sop)[$(x.op)]")
