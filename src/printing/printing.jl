struct TextStyle end
struct LatexStyle end
const PrintStyle = Union{TextStyle, LatexStyle}

function Base.show(io::IO, h::Union{FockSpace, NLevelSpace, CollectiveNLevelSpace, PauliSpace, SpinSpace, PhaseSpace})
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

Base.show(io::IO, idx::Index) = print(io, index_name(idx))

function show_index_suffix(io::IO, idx::Index)
    if has_index(idx)
        write(io, "_")
        print(io, index_name(idx))
    end
    return
end

const subscript_digits = ('₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉')
function write_subscript(io::IO, n::Int)
    return if 0 <= n <= 9
        write(io, subscript_digits[n + 1])
    else
        for d in reverse(digits(n))
            write(io, subscript_digits[d + 1])
        end
    end
end

const xyz_sym = (:x, :y, :z)

function Base.show(io::IO, x::Op)
    print(io, operator_name(x))
    show_index_suffix(io, x.index)
    k = x.kind
    if k === OP_CREATE
        write(io, "'")
    elseif k === OP_TRANSITION || k === OP_COLLECTIVE_TRANSITION
        write_subscript(io, Int(x.l1))
        write_subscript(io, Int(x.l2))
    elseif k === OP_PAULI || k === OP_SPIN
        print(io, xyz_sym[x.l1])
    end
    return
end

is_call_to(x, heads::Tuple) = SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) in heads

needs_im_parens(v::Num) = SymbolicUtils.iscall(SymbolicUtils.unwrap(v))
is_loose_head(v::Num) = is_call_to(SymbolicUtils.unwrap(v), (/, +))

function show_part(io::IO, v, brace::Bool)
    brace && write(io, "(")
    print(io, v)
    brace && write(io, ")")
    return
end

function show_prefactor(io::IO, c::CNum)
    tail = c.tail
    tail isa RawSymbolicCoeff && return print(io, tail.expr)
    return show_display(io, to_num(c))
end

function show_display(io::IO, c::Complex{Num})
    return if iszero_num(imag(c))
        print(io, real(c))
    elseif iszero_num(real(c))
        i = imag(c)
        if isone(i)
            write(io, "im")
        elseif isequal(i, Num(-1))
            write(io, "-im")
        else
            show_part(io, i, needs_im_parens(i))
            write(io, "im")
        end
    else
        re, i = real(c), imag(c)
        write(io, "(")
        show_part(io, re, is_loose_head(re))
        write(io, " + ")
        show_part(io, i, needs_im_parens(i))
        write(io, "im)")
    end
end

is_unit(c::CNum) = isequal(c, CNUM_ONE)
is_neg_unit(c::CNum) = isequal(c, CNUM_NEG1)

function is_real_negative(c::CNum)
    is_native(c) && return imag(native_scalar(c)) == 0 && real(native_scalar(c)) < 0
    t = c.tail
    t isa Poly &&
        return all_real_negative(t.terms)
    return is_real_negative_sym(c)
end
all_real_negative(terms::Vector{<:Monomial}) = all(terms) do m
    scalar = monomial_scalar(m)
    imag(scalar) == 0 && real(scalar) < 0
end
@noinline function is_real_negative_sym(c::CNum)
    re, im = realimag(c)
    iszero_num(im) || return false
    return const_sign(SymbolicUtils.unwrap(re)) < 0
end
needs_pf_parens(c::Complex{Num}) = iszero(imag(c)) && is_loose_head(real(c))

function show_scaled_prefactor(io::IO, c::CNum)
    is_neg_unit(c) && return write(io, "-")
    is_unit(c) && return
    tail = c.tail
    if tail isa RawSymbolicCoeff
        show_part(io, tail.expr, is_call_to(tail.expr, (+, /)))
    else
        d = to_num(c)
        brace = !is_native(c) && needs_pf_parens(d)
        brace && write(io, "(")
        show_display(io, d)
        brace && write(io, ")")
    end
    write(io, " * ")
    return
end

function show_ops(io::IO, ops::Vector{Op})
    for (k, op) in enumerate(ops)
        k > 1 && write(io, " * ")
        show(io, op)
    end
    return
end

function show_term(io::IO, c::CNum, ops::Vector{Op})
    isempty(ops) && return show_prefactor(io, c)
    show_scaled_prefactor(io, c)
    show_ops(io, ops)
    return
end

function term_signature(t::QAdd)
    term, c = first(t.arguments)
    return term.ops, c, term.ne
end

function show_terms(io::IO, st::Vector{QAdd})
    isempty(st) && return write(io, "0")
    ops1, c1, _ = term_signature(st[1])
    show_term(io, c1, ops1)
    for i in 2:length(st)
        ops_i, c_i, _ = term_signature(st[i])
        if is_real_negative(c_i)
            write(io, " - ")
            show_term(io, -c_i, ops_i)
        else
            write(io, " + ")
            show_term(io, c_i, ops_i)
        end
    end
    return
end

function group_dep_terms(dep_qadds::Vector{QAdd}, indices::Vector{Index})
    groups = Tuple{Vector{Index}, Vector{NonEqualPair}, Vector{QAdd}}[]
    for q in dep_qadds
        ops, c, ne = term_signature(q)
        used = bound_indices(c, ops, ne, indices)
        slot = position_of(g -> isequal(g[1], used) && isequal(g[2], ne), groups)
        if iszero(slot)
            push!(groups, (used, copy_ne(ne), QAdd[q]))
        else
            push!(groups[slot][3], q)
        end
    end
    return groups
end

function show_sum_prefix(io::IO, ::TextStyle, indices::Vector{Index}, ne_pairs::Vector{NonEqualPair})
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
    write(io, " ")
    return
end

function show_sum_group(io::IO, terms::Vector{QAdd}, indices::Vector{Index}, ne_pairs::Vector{NonEqualPair})
    show_sum_prefix(io, TextStyle(), indices, ne_pairs)
    group = length(terms) > 1
    group && write(io, "(")
    show_terms(io, terms)
    group && write(io, ")")
    return io
end

function SymbolicUtils.show_call(io::IO, ::typeof(expim), x::SymbolicUtils.BasicSymbolic; kw...)
    neg, body = split_phase_sign(only(SymbolicUtils.arguments(x)))
    write(io, neg ? "exp(-im*" : "exp(im*")
    show_part(io, body, is_call_to(body, (+,)))
    write(io, ")")
    return io
end

function split_phase_sign(arg)
    neg = leading_sign(arg) < 0
    return neg, neg ? SymbolicUtils.unwrap(expand(-arg)) : arg
end

function SymbolicUtils.show_call(io::IO, ::SumFunc, x::SymbolicUtils.BasicSymbolic; kw...)
    show_sum_prefix(io, TextStyle(), sum_indices(x), sum_ne(x))
    body = sum_body(x)
    show_part(io, body, is_call_to(body, (+,)))
    return io
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
    if isempty(x.indices)
        show_terms(io, st)
        return
    end
    dep, indep = split_index_dependent(st, x.indices)
    isempty(indep) || show_terms(io, indep)
    isempty(indep) || isempty(dep) || write(io, " + ")
    isempty(dep) || show_sum_groups(io, dep, x.indices)
    return
end

function split_index_dependent(st::Vector{QAdd}, indices::Vector{Index})
    dep = QAdd[]
    indep = QAdd[]
    for t in st
        ops, c, ne = term_signature(t)
        dependent = any(idx -> depends_on_scope_index(c, ops, ne, idx), indices)
        push!(dependent ? dep : indep, t)
    end
    return dep, indep
end

function show_sum_groups(io::IO, dep::Vector{QAdd}, indices::Vector{Index})
    for (k, (used, ne_pairs, terms)) in enumerate(group_dep_terms(dep, indices))
        k > 1 && write(io, " + ")
        show_sum_group(io, terms, used, ne_pairs)
    end
    return
end

const SHOW_RULE_LIMIT = 4

function show_unitary_rule(io::IO, U::UnitaryTransform, generator::Op)
    show(io, generator)
    write(io, " ↦ ")
    return show(io, U.rules[generator])
end

rules_elided(io::IO, U::UnitaryTransform) = get(io, :limit, false) && length(U.rules) > SHOW_RULE_LIMIT

function Base.show(io::IO, U::UnitaryTransform)
    write(io, "UnitaryTransform(")
    if rules_elided(io, U)
        print(io, length(U.rules), " rules")
    else
        for (k, generator) in enumerate(generators(U))
            k > 1 && write(io, ", ")
            show_unitary_rule(io, U, generator)
        end
    end
    timed = U isa UnitaryTransform{DynamicTime}
    timed && print(io, "; time = ", U.time.variable)
    if !iszero(U.gauge)
        write(io, timed ? ", gauge = " : "; gauge = ")
        show(io, U.gauge)
    end
    return write(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", U::UnitaryTransform)
    n = length(U.rules)
    if U isa UnitaryTransform{DynamicTime}
        write(io, "Time-dependent UnitaryTransform in ")
        show(io, U.time.variable)
    else
        write(io, "UnitaryTransform")
    end
    print(io, " with ", n, n == 1 ? " rule" : " rules")
    if !rules_elided(io, U)
        write(io, ":")
        for generator in generators(U)
            write(io, "\n  ")
            show_unitary_rule(io, U, generator)
        end
    end
    if !iszero(U.gauge)
        write(io, "\n\nGauge:\n  ")
        show(io, U.gauge)
    end
    return io
end
