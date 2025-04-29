using Symbolics
import Symbolics: simplify,Term
using SymbolicUtils
import SymbolicUtils: Symbolic,_isone,flatten_term,isnotflat,Chain,Fixpoint,Prewalk,sorted_arguments
using TermInterface
import TermInterface: isexpr,head,iscall,children,operation,arguments,metadata,maketerm
using Moshi.Data: @data, variant_type, variant_name
using Moshi.Match: @match

using LinearAlgebra
import LinearAlgebra: eigvecs,ishermitian,conj,transpose,inv,exp,vec,tr

import QuantumInterface:
    apply!,
    tensor, ⊗,
    basis,Basis,samebases,IncompatibleBases,SpinBasis,FockBasis,CompositeBasis,
    nqubits,
    projector,dagger,tr,ptrace,
    AbstractBra,AbstractKet,AbstractOperator,AbstractSuperOperator,
    express,AbstractRepresentation,AbstractUse,UseAsState,UseAsObservable,UseAsOperation,
    QuantumOpticsRepr,QuantumMCRepr,CliffordRepr

export SymQObj,QObj,
       AbstractRepresentation,AbstractUse,
       QuantumOpticsRepr,QuantumMCRepr,CliffordRepr,
       UseAsState,UseAsObservable,UseAsOperation,
       apply!,
       express,
       tensor,⊗,
       dagger,projector,commutator,anticommutator,conj,transpose,inv,exp,vec,tr,ptrace,
       I,X,Y,Z,σˣ,σʸ,σᶻ,Pm,Pp,σ₋,σ₊,
       H,CNOT,CPHASE,XCX,XCY,XCZ,YCX,YCY,YCZ,ZCX,ZCY,ZCZ,
       X1,X2,Y1,Y2,Z1,Z2,X₁,X₂,Y₁,Y₂,Z₁,Z₂,L0,L1,Lp,Lm,Lpi,Lmi,L₀,L₁,L₊,L₋,L₊ᵢ,L₋ᵢ,
       vac,F₀,F0,F₁,F1,inf_fock_basis,
       N,n̂,Create,âꜛ,Destroy,â,basis,SpinBasis,FockBasis,
       SBra,SKet,SOperator,SHermitianOperator,SUnitaryOperator,SHermitianUnitaryOperator,SSuperOperator,
       @ket,@bra,@op,@superop,
       SAdd,SAddBra,SAddKet,SAddOperator,
       SScaled,SScaledBra,SScaledOperator,SScaledKet,
       STensorBra,STensorKet,STensorOperator,
       SZeroBra,SZeroKet,SZeroOperator,
       SConjugate,STranspose,SProjector,SDagger,SInvOperator,SExpOperator,SVec,STrace,SPartialTrace,
       MixedState,IdentityOp,
       SApplyKet,SApplyBra,SMulOperator,SSuperOpApply,SCommutator,SAnticommutator,SBraKet,SOuterKetBra,
       HGate,XGate,YGate,ZGate,CPHASEGate,CNOTGate,
       XBasisState,YBasisState,ZBasisState,FockState,CoherentState,SqueezedState,
       NumberOp,CreateOp,DestroyOp,PhaseShiftOp,DisplaceOp,SqueezeOp,
       XCXGate,XCYGate,XCZGate,YCXGate,YCYGate,YCZGate,ZCXGate,ZCYGate,ZCZGate,
       qsimplify,qsimplify_pauli,qsimplify_commutator,qsimplify_anticommutator,qsimplify_fock,
       qexpand,
       isunitary,
       KrausRepr,kraus

include("types.jl")

##
# Metadata cache helpers
##

const CacheType = Dict{Tuple{<:AbstractRepresentation,<:AbstractUse},Any}
mutable struct Metadata
    express_cache::CacheType # TODO use more efficient mapping
end
Metadata() = Metadata(CacheType())

"""Decorate a struct definition in order to add a metadata dict which would be storing cached `express` results."""
macro withmetadata(strct)
    ex = quote $strct end
    if @capture(ex, (struct T_{params__} fields__ end) | (struct T_{params__} <: A_ fields__ end))
        struct_name = namify(T)
        args = (namify(i) for i in fields if !MacroTools.isexpr(i, String, :string))
        constructor = :($struct_name{S}($(args...)) where S = new{S}($((args..., :(Metadata()))...)))
    elseif @capture(ex, struct T_ fields__ end)
        struct_name = namify(T)
        args = (namify(i) for i in fields if !MacroTools.isexpr(i, String, :string))
        constructor = :($struct_name($(args...)) = new($((args..., :(Metadata()))...)))
    else @capture(ex, struct T_ end)
        struct_name = namify(T)
        constructor = :($struct_name() = new($:(Metadata())))
    end
    struct_args = strct.args[end].args
    push!(struct_args, constructor, :(metadata::Metadata))
    esc(quote
    Base.@__doc__ $strct
    metadata(x::$struct_name)=x.metadata
    end)
end

##
# Utilities
##

include("utils.jl")

##
# Most symbolic objects defined here
##

include("literal_objects.jl")
include("basic_ops_homogeneous.jl")
include("basic_ops_inhomogeneous.jl")
include("basic_superops.jl")
include("linalg.jl")
include("predefined.jl")
include("predefined_CPTP.jl")
include("predefined_fock.jl")

##
# Symbolic and simplification rules
##

include("rules.jl")

##
# Expressing in specific formalism
##

include("express.jl")

##
# Printing
##

include("latexify.jl")
