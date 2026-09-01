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

function _show_prefactor(io::IO, c::CNum)
    tail = c.tail
    tail isa RawSymbolicCoeff && return print(io, tail.expr)
    return _show_display(io, to_num(c))
end

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
    re, im = _realimag(c)
    _iszero_num(im) || return false
    r = Symbolics.value(SymbolicUtils.unwrap(re))
    return r isa Real && r < 0
end
_is_real_negative(::Number) = false

# Only `_show_prefactor`'s real-only branch can leave a loose head exposed at top level;
# the pure-imaginary and mixed branches brace their own parts.
_needs_pf_parens(c::Complex{Num}) = iszero(imag(c)) && _is_loose_head(real(c))

function _show_term(io::IO, c::CNum, ops::Vector{Op})
    if isempty(ops)
        _show_prefactor(io, c)
        return
    end
    if _is_neg_unit(c)
        write(io, "-")
    elseif !_is_unit(c)
        tail = c.tail
        if tail isa RawSymbolicCoeff
            expr = tail.expr
            brace = SymbolicUtils.iscall(expr) &&
                (SymbolicUtils.operation(expr) === (+) || SymbolicUtils.operation(expr) === (/))
            brace && write(io, "(")
            print(io, expr)
            brace && write(io, ")")
        else
            # Lower once so the parenthesis decision and rendering inspect the same expression.
            d = to_num(c)
            brace = !_is_native(c) && _needs_pf_parens(d)
            brace && write(io, "(")
            _show_display(io, d)
            brace && write(io, ")")
        end
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

function _show_unitary_rule(io::IO, U::UnitaryTransform, generator::Op)
    show(io, generator)
    write(io, " ↦ ")
    return show(io, U.rules[generator])
end

_show_unitary_time(::IO, ::UnitaryTransform{StaticTime}) = false
function _show_unitary_time(io::IO, U::UnitaryTransform{DynamicTime})
    write(io, "time = ")
    show(io, U.time.variable)
    return true
end

function _show_unitary_metadata(io::IO, U::UnitaryTransform)
    shown = _show_unitary_time(io, U)
    if !iszero(U.gauge)
        shown && write(io, ", ")
        write(io, "gauge = ")
        show(io, U.gauge)
        shown = true
    end
    return shown
end

function Base.show(io::IO, U::UnitaryTransform)
    gs = U.generators
    n = length(gs)
    write(io, "UnitaryTransform(")
    limited = get(io, :limit, false) && n > _SHOW_RULE_LIMIT
    if limited
        print(io, n, " rules")
    else
        for (k, generator) in enumerate(gs)
            k > 1 && write(io, ", ")
            _show_unitary_rule(io, U, generator)
        end
    end
    if U isa UnitaryTransform{DynamicTime} || !iszero(U.gauge)
        write(io, "; ")
        _show_unitary_metadata(io, U)
    end
    return write(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", U::UnitaryTransform)
    gs = U.generators
    n = length(gs)
    if U isa UnitaryTransform{DynamicTime}
        write(io, "Time-dependent UnitaryTransform in ")
        show(io, U.time.variable)
    else
        write(io, "UnitaryTransform")
    end
    print(io, " with ", n, n == 1 ? " rule" : " rules")

    limited = get(io, :limit, false) && n > _SHOW_RULE_LIMIT
    if !limited
        write(io, ":")
        for generator in gs
            write(io, "\n  ")
            _show_unitary_rule(io, U, generator)
        end
    end
    if !iszero(U.gauge)
        write(io, "\n\nGauge:\n  ")
        show(io, U.gauge)
    end
    return nothing
end
