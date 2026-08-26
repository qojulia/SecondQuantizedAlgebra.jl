module ClosurePrototype

using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: Coeff, HilbertSpace, Op, QAdd, QField, QSym, QTerm,
    _CNUM_ZERO

"""Deterministic work limits for disposable closure discovery."""
Base.@kwdef struct ClosureLimits
    max_basis::Int = 128
    max_degree::Int = 8
    max_commutators::Int = 256
end

"""
Transient ordered basis. `QTerm` is reused deliberately: it already represents an exact
operator monomial together with its inequality constraints and cached structural hash.
"""
struct AdjointBasisPrototype{H <: HilbertSpace}
    space::H
    terms::Vector{QTerm}
    lookup::Dict{QTerm, Int32}
    max_degree::Int
end

struct ClosureResult{H <: HilbertSpace}
    basis::AdjointBasisPrototype{H}
    action::Matrix{Coeff}
    offset::Vector{Coeff}
    commutators::Int
    status::Symbol
    offending_degree::Int
end

"""Competing wrapper used only to quantify whether duplicating `QTerm` buys anything."""
struct BasisMonomial
    ops::Vector{Op}
    ne::Vector{Tuple{SecondQuantizedAlgebra.Index, SecondQuantizedAlgebra.Index}}
    hash::UInt
end

@enum InferredFamily::UInt8 begin
    FAMILY_UNKNOWN
    FAMILY_FOCK
    FAMILY_PHASE
    FAMILY_PAULI
    FAMILY_SPIN
    FAMILY_NLEVEL
    FAMILY_COLLECTIVE_NLEVEL
    FAMILY_CONFLICT
end

"""Concrete metadata available from one occupied `space_index`."""
struct InferredSlot
    family::InferredFamily
    complete::Bool
    n::Int32
    ground::Int32
    names::Vector{Symbol}
    indexed::Bool
end

struct InferredContext
    slots::Dict{Int32, InferredSlot}
    contiguous::Bool
end

BasisMonomial(t::QTerm) = BasisMonomial(copy(t.ops), copy(t.ne), hash(t))
Base.hash(t::BasisMonomial, h::UInt) = hash(t.hash, h)
Base.isequal(a::BasisMonomial, b::BasisMonomial) =
    a.hash == b.hash && isequal(a.ops, b.ops) && isequal(a.ne, b.ne)
Base.:(==)(a::BasisMonomial, b::BasisMonomial) = isequal(a, b)

_as_qadd(x::QAdd) = x
_as_qadd(x::QSym) = 1 * x

function _term_expression(term::QTerm)
    isempty(term.ne) || throw(ArgumentError("prototype cannot reconstruct scoped inequalities"))
    isempty(term.ops) && return 1
    return foldl(*, term.ops)
end

function _normalized_commutator(generator::QAdd, term::QTerm)
    q = expand_completeness(commutator(generator, _term_expression(term)))
    isempty(q.indices) || throw(ArgumentError("summation-scoped closure requires a separate indexed prototype"))
    return q
end

function _seed_terms(seeds)
    result = QTerm[]
    lookup = Dict{QTerm, Int32}()
    for seed in seeds
        for candidate in (seed, adjoint(seed))
            q = expand_completeness(_as_qadd(candidate))
            isempty(q.indices) || throw(ArgumentError("summation-scoped seed is unsupported"))
            for (term, _) in q
                isempty(term.ops) && continue
                if !haskey(lookup, term)
                    push!(result, term)
                    lookup[term] = Int32(length(result))
                end
            end
        end
    end
    return result, lookup
end

function discover_closure(
        space::H, generator, seeds;
        limits::ClosureLimits = ClosureLimits(),
    ) where {H <: HilbertSpace}
    G = _as_qadd(generator)
    terms, lookup = _seed_terms(seeds)
    maxdegree = isempty(terms) ? 0 : maximum(length(t.ops) for t in terms)
    columns = QAdd[]
    status = :closed
    offending_degree = 0
    cursor = 1
    commutators = 0

    while cursor <= length(terms)
        if commutators >= limits.max_commutators
            status = :max_commutators
            break
        end
        q = try
            _normalized_commutator(G, terms[cursor])
        catch error
            error isa ArgumentError || rethrow()
            status = :indexed_scope
            break
        end
        commutators += 1
        push!(columns, q)
        for (term, _) in q
            isempty(term.ops) && continue
            degree = length(term.ops)
            maxdegree = max(maxdegree, degree)
            if degree > limits.max_degree
                status = :max_degree
                offending_degree = degree
                break
            end
            if !haskey(lookup, term)
                if length(terms) >= limits.max_basis
                    status = :max_basis
                    offending_degree = degree
                    break
                end
                push!(terms, term)
                lookup[term] = Int32(length(terms))
            end
        end
        status === :closed || break
        cursor += 1
    end

    n = length(terms)
    action = fill(_CNUM_ZERO, n, n)
    offset = fill(_CNUM_ZERO, n)
    for column in eachindex(columns)
        for (term, coefficient) in columns[column]
            if isempty(term.ops)
                offset[column] = coefficient
            else
                row = get(lookup, term, Int32(0))
                row == 0 && continue
                action[row, column] = coefficient
            end
        end
    end
    basis = AdjointBasisPrototype(space, terms, lookup, maxdegree)
    return ClosureResult(basis, action, offset, commutators, status, offending_degree)
end

"""Infer only algebra metadata present in `Op`; this cannot recreate a Hilbert space."""
function infer_operator_context(operators::AbstractVector{Op})
    slots = Dict{Int32, InferredSlot}()
    for op in operators
        family, complete, n, ground = if is_destroy(op) || is_create(op)
            (FAMILY_FOCK, true, Int32(0), Int32(0))
        elseif is_position(op) || is_momentum(op)
            (FAMILY_PHASE, true, Int32(0), Int32(0))
        elseif is_pauli(op)
            (FAMILY_PAULI, true, Int32(0), Int32(0))
        elseif is_spin(op)
            (FAMILY_SPIN, true, Int32(0), Int32(0))
        elseif is_transition(op)
            (FAMILY_NLEVEL, true, op.nlev, op.g)
        elseif is_collective_transition(op)
            (FAMILY_COLLECTIVE_NLEVEL, false, Int32(0), Int32(0))
        else
            (FAMILY_UNKNOWN, false, Int32(0), Int32(0))
        end
        item = InferredSlot(
            family, complete, n, ground, Symbol[operator_name(op)],
            op.index !== SecondQuantizedAlgebra.NO_INDEX,
        )
        existing = get(slots, op.space_index, nothing)
        if existing === nothing
            slots[op.space_index] = item
        elseif existing.family != family
            slots[op.space_index] = InferredSlot(
                FAMILY_CONFLICT, false, 0, 0,
                unique(vcat(existing.names, item.names)),
                existing.indexed || item.indexed,
            )
        else
            slots[op.space_index] = InferredSlot(
                existing.family, existing.complete,
                max(existing.n, item.n), max(existing.ground, item.ground),
                unique(vcat(existing.names, item.names)),
                existing.indexed || item.indexed,
            )
        end
    end
    slot_ids = sort!(collect(keys(slots)))
    contiguous = isempty(slot_ids) || slot_ids == collect(Int32(1):maximum(slot_ids))
    return InferredContext(slots, contiguous)
end

function warm_measure(f; samples::Int = 21)
    f()
    times = Float64[]
    allocs = Int[]
    sizehint!(times, samples)
    sizehint!(allocs, samples)
    for _ in 1:samples
        GC.gc()
        push!(allocs, @allocated f())
        push!(times, @elapsed f())
    end
    sort!(times)
    sort!(allocs)
    middle = (samples + 1) >> 1
    return (; seconds = times[middle], bytes = allocs[middle])
end

function storage_comparison(terms::Vector{QTerm}; repetitions::Int = 1000)
    qdict = Dict(term => Int32(i) for (i, term) in pairs(terms))
    wrappers = BasisMonomial.(terms)
    bdict = Dict(term => Int32(i) for (i, term) in pairs(wrappers))
    qconstruct = warm_measure(() -> begin
        for _ in 1:repetitions
            Dict(term => Int32(i) for (i, term) in pairs(terms))
        end
    end; samples = 11)
    bconstruct = warm_measure(() -> begin
        for _ in 1:repetitions
            converted = BasisMonomial.(terms)
            Dict(term => Int32(i) for (i, term) in pairs(converted))
        end
    end; samples = 11)
    qlookup = warm_measure(() -> begin
        total = 0
        for _ in 1:repetitions, term in terms
            total += qdict[term]
        end
        total
    end; samples = 11)
    blookup = warm_measure(() -> begin
        total = 0
        for _ in 1:repetitions, term in wrappers
            total += bdict[term]
        end
        total
    end; samples = 11)
    return (; qconstruct, bconstruct, qlookup, blookup)
end

end

