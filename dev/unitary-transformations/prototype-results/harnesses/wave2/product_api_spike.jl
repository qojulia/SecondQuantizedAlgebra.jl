module ProductAPIPrototype

using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: CNum, QAdd, QTerm, Op, HilbertSpace, ProductSpace,
    operator_name, operator_index, optype, index_name, index_range, index_slot,
    is_destroy, is_create, is_transition, is_collective_transition, is_pauli,
    is_spin, is_position, is_momentum

export FlatCoordinates, AnalysisResult, analyze_frame, reset_analysis_count!,
    analysis_count, project_action, alpha_signature, validate_hilbert

"""Disposable flat coordinate store; support is analysis metadata, not key identity."""
struct FlatCoordinates{H <: HilbertSpace}
    hilbert::H
    terms::Vector{QTerm}
    lookup::Dict{QTerm, Int32}
    support::Vector{Vector{Int32}}
end

struct AnalysisResult{H <: HilbertSpace}
    hilbert::H
    coordinates::FlatCoordinates{H}
end

const _ANALYSIS_COUNT = Ref(0)
reset_analysis_count!() = (_ANALYSIS_COUNT[] = 0)
analysis_count() = _ANALYSIS_COUNT[]

@inline _factor(h::ProductSpace, si::Int) = h.spaces[si]
@inline _factor(h::HilbertSpace, si::Int) = si == 1 ? h : throw(ArgumentError("operator targets factor $si in a one-factor Hilbert space"))

function _role_matches(op::Op, factor::HilbertSpace)
    (is_destroy(op) || is_create(op)) && return factor isa FockSpace
    is_transition(op) && return factor isa NLevelSpace
    is_collective_transition(op) && return factor isa CollectiveNLevelSpace
    is_pauli(op) && return factor isa PauliSpace
    is_spin(op) && return factor isa SpinSpace
    (is_position(op) || is_momentum(op)) && return factor isa PhaseSpace
    return false
end

function _validate_op(h::HilbertSpace, op::Op)
    si = Int(op.space_index)
    1 <= si <= length(h) || throw(ArgumentError("operator $(operator_name(op)) targets absent factor $si"))
    factor = _factor(h, si)
    _role_matches(op, factor) || throw(ArgumentError("operator $(operator_name(op)) is incompatible with factor $si ($(typeof(factor)))"))
    idx = operator_index(op)
    if has_index(idx)
        Int(idx.space_index) == si || throw(ArgumentError("operator/index factor mismatch: operator $si, index $(idx.space_index)"))
        index_slot(idx) === nothing &&
            throw(ArgumentError("free indexed families are not supported by the initial exact-closure boundary"))
    end
    if is_transition(op)
        (Int(op.nlev) == factor.n && Int(op.g) == factor.ground_state) ||
            throw(ArgumentError("transition metadata disagrees with factor $si"))
    elseif is_collective_transition(op)
        max(Int(op.l1), Int(op.l2)) <= factor.n ||
            throw(ArgumentError("collective-transition level exceeds factor $si"))
    end
    return nothing
end

function validate_hilbert(h::HilbertSpace, expressions...)
    for expr in expressions
        expr isa Op && (_validate_op(h, expr); continue)
        expr isa QAdd || throw(ArgumentError("expected operator expression, got $(typeof(expr))"))
        for (term, _) in expr
            for op in term.ops
                _validate_op(h, op)
            end
        end
    end
    return h
end

function _terms(expr::QAdd)
    out = QTerm[]
    sizehint!(out, length(expr))
    for (term, _) in expr
        isempty(term.ops) || push!(out, term)
    end
    return out
end
_terms(op::Op) = _terms(op + 0)

function _support(term::QTerm)
    s = Int32[]
    for op in term.ops
        si = Int32(op.space_index)
        si in s || push!(s, si)
    end
    sort!(s)
    return s
end

function FlatCoordinates(h::H, expressions...) where {H <: HilbertSpace}
    validate_hilbert(h, expressions...)
    terms = QTerm[]
    lookup = Dict{QTerm, Int32}()
    support = Vector{Vector{Int32}}()
    for expr in expressions, term in _terms(expr)
        haskey(lookup, term) && continue
        push!(terms, term)
        lookup[term] = Int32(length(terms))
        push!(support, _support(term))
    end
    return FlatCoordinates{H}(h, terms, lookup, support)
end

"""Project an exact commutator into flat QTerm coordinates, refusing growth."""
function project_action(G, coordinates::FlatCoordinates)
    n = length(coordinates.terms)
    action = Matrix{CNum}(undef, n, n)
    fill!(action, 0)
    for j in 1:n
        X = prod(coordinates.terms[j].ops)
        Y = im * commutator(G, X)
        for (term, coefficient) in Y
            isempty(term.ops) && continue
            row = get(coordinates.lookup, term, Int32(0))
            row == 0 && throw(ArgumentError("closure escaped at coordinate $j through $term"))
            action[Int(row), j] = coefficient
        end
    end
    return action
end

# Thin public-shape adapter. Keyword handling ends before private positional analysis.
"""
    analyze_frame(generator::QAdd; hilbert)

Disposable public-shape experiment: validate the explicit Hilbert space, then dispatch to a
concrete private positional method. Exact Hilbert identity is not reconstructible from `Op`.
"""
function analyze_frame(generator::QAdd; hilbert = nothing)
    hilbert isa HilbertSpace ||
        throw(ArgumentError("`hilbert` is required when the complete Hilbert space cannot be inferred"))
    validate_hilbert(hilbert, generator) # mismatch is rejected before the counter/analysis
    return _analyze_frame(generator, hilbert)
end

_adapter_call(generator::QAdd, hilbert::H) where {H <: HilbertSpace} =
    analyze_frame(generator; hilbert)

function _analyze_frame(generator, hilbert::H) where {H <: HilbertSpace}
    _ANALYSIS_COUNT[] += 1
    coordinates = FlatCoordinates(hilbert, generator)
    return AnalysisResult{H}(hilbert, coordinates)
end

# A future family descriptor can canonicalize bound labels before exact Dict lookup.
# This is deliberately only a signature experiment, not wildcard QTerm equality.
function alpha_signature(term::QTerm)
    binders = Dict{Tuple{Int32, Int32, Int32}, Int32}()
    next = Int32(0)
    ops = map(term.ops) do op
        idx = operator_index(op)
        binder = Int32(0)
        if has_index(idx)
            key = (idx.name_id, idx.range_id, idx.space_index)
            binder = get!(binders, key) do
                next += 1
                next
            end
        end
        (optype(op), op.space_index, binder, op.l1, op.l2, op.g, op.nlev)
    end
    inequalities = Tuple{Int32, Int32}[]
    for (a, b) in term.ne
        ka = (a.name_id, a.range_id, a.space_index)
        kb = (b.name_id, b.range_id, b.space_index)
        haskey(binders, ka) && haskey(binders, kb) || continue
        ia, ib = binders[ka], binders[kb]
        push!(inequalities, ia < ib ? (ia, ib) : (ib, ia))
    end
    sort!(unique!(inequalities))
    return (Tuple(ops), Tuple(inequalities), next)
end

end

