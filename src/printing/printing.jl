function Base.show(io::IO, h::FockSpace)
    write(io, "ℋ(")
    print(io, h.name)
    return write(io, ")")
end
function Base.show(io::IO, h::ProductSpace)
    show(io, h.spaces[1])
    for i in 2:length(h.spaces)
        write(io, " ⊗ ")
        show(io, h.spaces[i])
    end
    return
end

Base.show(io::IO, h::NLevelSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::CollectiveNLevelSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::PauliSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::SpinSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::PhaseSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))

Base.show(io::IO, idx::Index) = print(io, index_name(idx))

function _show_index_suffix(io::IO, idx::Index)
    if has_index(idx)
        write(io, "_")
        print(io, index_name(idx))
    end
    return
end

const _subscript_digits = ('₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉')
function _write_subscript(io::IO, n::Int)
    return if 0 <= n <= 9
        write(io, _subscript_digits[n + 1])
    else
        for d in reverse(digits(n))
            write(io, _subscript_digits[d + 1])
        end
    end
end

const _xyz_sym = (:x, :y, :z)

function Base.show(io::IO, x::Op)
    print(io, operator_name(x))
    _show_index_suffix(io, x.index)
    k = x.kind
    if k === OP_CREATE
        write(io, "'")
    elseif k === OP_TRANSITION || k === OP_COLLECTIVE_TRANSITION
        _write_subscript(io, Int(x.l1))
        _write_subscript(io, Int(x.l2))
    elseif k === OP_PAULI || k === OP_SPIN
        print(io, _xyz_sym[x.l1])
    end
    return
end

# The `im` suffix is juxtaposition, so every call has to be braced, not just the loose heads:
# `x^2` reads as `x^(2im)`, `x*y` merges into the identifier `yim`, `sqrt(x)im` is a syntax
# error.
function _needs_im_parens(v::Num)
    u = SymbolicUtils.unwrap(v)
    return u isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(u)
end

# A `/` or `+` head has looser precedence than the ` * ops` that follows a prefactor, so
# such a part has to be parenthesized.
function _is_loose_head(v::Num)
    _needs_im_parens(v) || return false
    op = SymbolicUtils.operation(SymbolicUtils.unwrap(v))
    return op === (/) || op === (+)
end

function _show_part(io::IO, v::Num, brace::Bool)
    brace && write(io, "(")
    print(io, v)
    brace && write(io, ")")
    return
end

# A display-pair key deliberately excludes the phase sign and scalar. Two monomials match
# when they contain the same phase atom at opposite nonzero integer powers and otherwise
# have identical factors. The lookup is display-only; none of these lowered forms are ever
# fed back into coefficient recognition.
struct _PhaseDisplayKey
    phase_id::UInt
    power::Int
    factor_ids::Vector{UInt}
    exps::Vector{Rational{Int}}
end

function Base.isequal(a::_PhaseDisplayKey, b::_PhaseDisplayKey)
    a.phase_id == b.phase_id || return false
    a.power == b.power || return false
    length(a.factor_ids) == length(b.factor_ids) || return false
    @inbounds for i in eachindex(a.factor_ids)
        a.factor_ids[i] == b.factor_ids[i] || return false
        a.exps[i] == b.exps[i] || return false
    end
    return true
end

function Base.hash(k::_PhaseDisplayKey, h::UInt)
    h = hash(k.power, hash(k.phase_id, h))
    @inbounds for i in eachindex(k.factor_ids)
        h = hash(k.exps[i], hash(k.factor_ids[i], h))
    end
    return h
end

struct _PhaseDisplayTerm{A, F}
    key::_PhaseDisplayKey
    positive::Bool
    argument::A
    factors::F
end

struct _PhasePending
    index::Int
    positive::Bool
end

# Describe a monomial with exactly one phase factor, or return `nothing`.
function _phase_display_term(m::Monomial)
    pos = 0
    @inbounds for i in eachindex(m.syms)
        if _is_phase(m.syms[i])
            pos != 0 && return nothing
            pos = i
        end
    end
    pos == 0 && return nothing
    e = m.exps[pos]
    denominator(e) == 1 || return nothing
    power = numerator(e)
    (iszero(power) || power == typemin(Int)) && return nothing
    argument = only(SymbolicUtils.arguments(m.syms[pos]))
    argument isa SymbolicUtils.BasicSymbolic || return nothing
    n = length(m.syms) - 1
    syms = Vector{SymbolicUtils.BasicSymbolic}(undef, n)
    factor_ids = Vector{UInt}(undef, n)
    exps = Vector{Rational{Int}}(undef, n)
    k = 0
    @inbounds for i in eachindex(m.syms)
        i == pos && continue
        k += 1
        syms[k] = m.syms[i]
        factor_ids[k] = objectid(m.syms[i])
        exps[k] = m.exps[i]
    end
    key = _PhaseDisplayKey(objectid(m.syms[pos]), abs(power), factor_ids, exps)
    return _PhaseDisplayTerm(key, power > 0, argument, syms)
end

function _phase_display_pairs(terms::Vector{Monomial})
    n = length(terms)
    partners = zeros(Int, n)
    descriptions = Vector{Union{Nothing, _PhaseDisplayTerm}}(undef, n)
    fill!(descriptions, nothing)
    pending = Dict{_PhaseDisplayKey, _PhasePending}()
    @inbounds for i in 1:n
        desc = _phase_display_term(terms[i])
        descriptions[i] = desc
        desc === nothing && continue
        previous = get(pending, desc.key, nothing)
        if previous === nothing
            pending[desc.key] = _PhasePending(i, desc.positive)
        elseif previous.positive != desc.positive
            partners[i] = previous.index
            partners[previous.index] = i
            delete!(pending, desc.key)
        end
    end
    return partners, descriptions
end

function _without_phase(desc::_PhaseDisplayTerm)
    return real(_term_to_num(Monomial(_ONE_C, desc.factors, desc.key.exps)))
end

# Display lowering. A conjugate pair `c₊p^n + c₋p^-n` is exactly
# `(c₊+c₋)cos(n*x) + im*(c₊-c₋)sin(n*x)`, so folding it keeps a rotation reading as a
# rotation instead of a sum of exponentials. Everything else lowers as usual, and nothing
# here feeds `_recognize`: the polynomial stays the representation.
function _display_num(p::Poly)
    n = length(p.terms)
    acc = Complex(_NUM_ZERO, _NUM_ZERO)
    partners, descriptions = _phase_display_pairs(p.terms)
    @inbounds for i in 1:n
        partner = partners[i]
        (partner != 0 && partner < i) && continue
        m = p.terms[i]
        if partner == 0
            acc += _term_to_num(m)
            continue
        end
        desc = descriptions[i]::_PhaseDisplayTerm
        cp = desc.positive ? m.scalar : p.terms[partner].scalar
        cm = desc.positive ? p.terms[partner].scalar : m.scalar
        rest = _without_phase(desc)
        θ = expand(Num(desc.key.power) * Num(desc.argument))
        # `cos` is even and `sin` odd, so a negative argument folds onto the sine's scalar
        # and `cos(-ω*t)` never reaches the page.
        d = cp - cm
        if _leading_sign(θ) < 0
            θ = expand(-θ)
            d = -d
        end
        acc += (
            Complex(_num_from_float(real(cp + cm)), _num_from_float(imag(cp + cm))) * cos(θ) +
                im * Complex(_num_from_float(real(d)), _num_from_float(imag(d))) * sin(θ)
        ) * rest
    end
    return acc
end

# `to_num` for display: identical except that phase pairs fold back to `cos`/`sin`.
_display_coeff(c::Coeff) =
    c.tail isa Poly ? Complex(_num_from_float(real(c.z)), _num_from_float(imag(c.z))) +
    _display_num(c.tail) : to_num(c)

_show_prefactor(io::IO, c::CNum) = _show_display(io, _display_coeff(c))

function _show_display(io::IO, c::Complex{Num})
    return if iszero(imag(c))
        # A loose head here is parenthesized by `_needs_pf_parens` at the call site,
        # which also lets a standalone constant term print without parentheses.
        print(io, real(c))
    elseif iszero(real(c))
        i = imag(c)
        if isone(i)
            write(io, "im")
        elseif isequal(i, Num(-1))
            write(io, "-im")
        else
            _show_part(io, i, _needs_im_parens(i))
            write(io, "im")
        end
    else
        re, i = real(c), imag(c)
        write(io, "(")
        _show_part(io, re, _is_loose_head(re))   # followed by ` + `, not by the suffix
        write(io, " + ")
        _show_part(io, i, _needs_im_parens(i))
        write(io, "im)")
    end
end
_show_prefactor(io::IO, c::Number) = print(io, c)

function _is_unit(c::CNum)
    return isequal(c, _CNUM_ONE)
end
_is_unit(c::Number) = isone(c)

function _is_neg_unit(c::CNum)
    return isequal(c, _CNUM_NEG1)
end
_is_neg_unit(c::Number) = isequal(c, -1)

function _is_real_negative(c::CNum)
    _is_native(c) && return imag(c.z) == 0 && real(c.z) < 0
    t = c.tail
    # A polynomial coefficient with every scalar real and negative factors out as ` - `,
    # which is most of what a rotation or a squeeze produces.
    t isa Poly && return all(m -> imag(m.scalar) == 0 && real(m.scalar) < 0, t.terms)
    return _is_real_negative_sym(c)
end
# Cold path: a non-native coefficient is a real negative only for symbolic constants
# that don't round-trip to `ComplexF64` (irrationals like `-π`, exact rationals).
# Isolated so the hot native branch stays concrete. `::Bool` pins the result: `<` on
# the abstract-typed `r` is unprovably `Bool` to inference (the Symbolics boundary),
# and leaving it `Any` would poison the `_show_terms` caller.
@noinline function _is_real_negative_sym(c::CNum)::Bool
    _iszero_num(imag(c)) || return false
    r = Symbolics.value(SymbolicUtils.unwrap(real(c)))
    return r isa Real && r < 0
end
_is_real_negative(::Number) = false

# Only `_show_prefactor`'s real-only branch can leave a loose head exposed at top level;
# the pure-imaginary and mixed branches brace their own parts.
# Judged on the displayed form, not the stored one: a folded phase pair prints as a bare
# `cos(...)` even though the polynomial behind it is a sum.
_needs_pf_parens(c::Complex{Num}) = iszero(imag(c)) && _is_loose_head(real(c))

function _show_term(io::IO, c::CNum, ops::Vector{Op})
    if isempty(ops)
        _show_prefactor(io, c)
        return
    end
    if _is_neg_unit(c)
        write(io, "-")
    elseif !_is_unit(c)
        # Lowered once: `_display_num`'s conjugate-pair pairing is O(n^2) and lowers every
        # monomial, so asking the parenthesis question and printing must share one result.
        d = _display_coeff(c)
        brace = !_is_native(c) && _needs_pf_parens(d)
        brace && write(io, "(")
        _show_display(io, d)
        brace && write(io, ")")
        write(io, " * ")
    end
    show(io, ops[1])
    for i in 2:length(ops)
        write(io, " * ")
        show(io, ops[i])
    end
    return
end

function _term_signature(t::QAdd)
    term, c = first(t.arguments)
    return term.ops, c, term.ne
end

function _show_terms(io::IO, st::Vector{QAdd})
    isempty(st) && return write(io, "0")
    ops1, c1, _ = _term_signature(st[1])
    _show_term(io, c1, ops1)
    for i in 2:length(st)
        ops_i, c_i, _ = _term_signature(st[i])
        if _is_real_negative(c_i)
            write(io, " - ")
            _show_term(io, -c_i, ops_i)
        else
            write(io, " + ")
            _show_term(io, c_i, ops_i)
        end
    end
    return
end

function _term_used_indices(t::QAdd, indices::Vector{Index})
    ops, c, _ = _term_signature(t)
    used = Index[]
    for idx in indices
        _depends_on_index_term(c, ops, idx) && push!(used, idx)
    end
    return used
end

function _group_dep_terms(dep_qadds, indices::Vector{Index})
    groups = Tuple{Vector{Index}, Vector{NonEqualPair}, Vector{QAdd}}[]
    for q in dep_qadds
        _, _, ne = _term_signature(q)
        used = _term_used_indices(q, indices)
        slot = findfirst(g -> isequal(g[1], used) && isequal(g[2], ne), groups)
        if slot === nothing
            push!(groups, (used, _copy_ne(ne), QAdd[q]))
        else
            push!(groups[slot][3], q)
        end
    end
    return groups
end

function _show_sum_prefix(io::IO, indices::Vector{Index}, ne_pairs::Vector{NonEqualPair})
    for idx in indices
        write(io, "Σ(")
        print(io, index_name(idx))
        write(io, "=1:")
        print(io, index_range(idx))
        write(io, ")")
    end
    if !isempty(ne_pairs)
        write(io, "(")
        for (k, (a, b)) in enumerate(ne_pairs)
            k > 1 && write(io, ",")
            print(io, index_name(a))
            write(io, "≠")
            print(io, index_name(b))
        end
        write(io, ")")
    end
    return nothing
end

function _show_sum_group(io::IO, terms::Vector{QAdd}, indices::Vector{Index}, ne_pairs::Vector{NonEqualPair})
    _show_sum_prefix(io, indices, ne_pairs)
    write(io, " ")
    if length(terms) > 1
        write(io, "(")
        _show_terms(io, terms)
        write(io, ")")
    else
        _show_terms(io, terms)
    end
    return nothing
end

# `expim(x)` reads as `exp(im*x)`, with a leading minus pulled onto the `im` so an inverse
# phase shows as `exp(-im*ω*t)` rather than `exp(im*(-ω*t))`.
function SymbolicUtils.show_call(io::IO, ::typeof(expim), x::SymbolicUtils.BasicSymbolic; kw...)
    arg = only(SymbolicUtils.arguments(x))
    neg = _leading_sign(arg) < 0
    write(io, neg ? "exp(-im*" : "exp(im*")
    # Negate on the `BasicSymbolic` itself: wrapping in `Num` first is a method error for a
    # non-`SymReal` vartype, which is what the argument of a phase can be.
    body = neg ? SymbolicUtils.unwrap(expand(-arg)) : arg
    paren = SymbolicUtils.iscall(body) && SymbolicUtils.operation(body) === (+)
    paren && write(io, "(")
    print(io, body)
    paren && write(io, ")")
    write(io, ")")
    return nothing
end

function SymbolicUtils.show_call(io::IO, ::SumFunc, x::SymbolicUtils.BasicSymbolic; kw...)
    _show_sum_prefix(io, _sum_indices(x), _sum_ne(x))
    write(io, " ")
    body = _sum_body(x)
    paren = SymbolicUtils.iscall(body) && SymbolicUtils.operation(body) === (+)
    paren && write(io, "(")
    print(io, body)
    paren && write(io, ")")
    return nothing
end

function SymbolicUtils.show_metadata(
        io::IO, x::SymbolicUtils.BasicSymbolic, ::Type{AverageOperator}, op,
    )
    write(io, "⟨")
    print(io, op)
    write(io, "⟩")
    if SymbolicUtils.iscall(x)
        write(io, "(")
        print(io, only(SymbolicUtils.arguments(x)))
        write(io, ")")
    end
    return true
end

function Base.show(io::IO, x::QAdd)
    st = sorted_arguments(x)
    isempty(st) && return write(io, "0")

    if !isempty(x.indices)
        dep = eltype(st)[]
        indep = eltype(st)[]
        for t in st
            ops, c, _ = _term_signature(t)
            if any(idx -> _depends_on_index_term(c, ops, idx), x.indices)
                push!(dep, t)
            else
                push!(indep, t)
            end
        end
        if !isempty(indep)
            _show_terms(io, indep)
            if !isempty(dep)
                write(io, " + ")
            end
        end
        if !isempty(dep)
            groups = _group_dep_terms(dep, x.indices)
            for (k, (used, ne_pairs, terms)) in enumerate(groups)
                k > 1 && write(io, " + ")
                _show_sum_group(io, terms, used, ne_pairs)
            end
        end
    else
        _show_terms(io, st)
    end
    return
end

# An N-level transform on 5 levels holds 25 rules of up to 25 terms each, so the full form is
# unreadable at a REPL prompt. `:limit` is what the REPL sets and a file write does not.
const _SHOW_RULE_LIMIT = 4

function Base.show(io::IO, U::UnitaryTransform)
    gs = generators(U)
    n = length(gs)
    shown = get(io, :limit, false) ? min(n, _SHOW_RULE_LIMIT) : n
    write(io, "UnitaryTransform(")
    for k in 1:shown
        k > 1 && write(io, ", ")
        show(io, gs[k])
        write(io, " → ")
        show(io, U.rules[gs[k]])
    end
    shown < n && write(io, ", … ($(n - shown) more)")
    if !iszero(U.gauge)
        write(io, "; gauge = ")
        show(io, U.gauge)
    end
    return write(io, ")")
end
