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

function needs_im_parens(v::Num)
    u = SymbolicUtils.unwrap(v)
    return u isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(u)
end

function is_loose_head(v::Num)
    needs_im_parens(v) || return false
    op = SymbolicUtils.operation(SymbolicUtils.unwrap(v))
    return op === (/) || op === (+)
end

function show_part(io::IO, v::Num, brace::Bool)
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
show_prefactor(io::IO, c::Number) = print(io, c)

function is_unit(c::CNum)
    return isequal(c, CNUM_ONE)
end
is_unit(c::Number) = isone(c)

function is_neg_unit(c::CNum)
    return isequal(c, CNUM_NEG1)
end
is_neg_unit(c::Number) = isequal(c, -1)

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
        expr = tail.expr
        brace = SymbolicUtils.iscall(expr) &&
            (SymbolicUtils.operation(expr) === (+) || SymbolicUtils.operation(expr) === (/))
        brace && write(io, "(")
        print(io, expr)
        brace && write(io, ")")
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

function show_term(io::IO, c::CNum, ops::Vector{Op})
    if isempty(ops)
        show_prefactor(io, c)
        return
    end
    show_scaled_prefactor(io, c)
    show(io, ops[1])
    for i in 2:length(ops)
        write(io, " * ")
        show(io, ops[i])
    end
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

function group_dep_terms(dep_qadds, indices::Vector{Index})
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

function show_sum_prefix(io::IO, indices::Vector{Index}, ne_pairs::Vector{NonEqualPair})
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
    return io
end

function show_sum_group(io::IO, terms::Vector{QAdd}, indices::Vector{Index}, ne_pairs::Vector{NonEqualPair})
    show_sum_prefix(io, indices, ne_pairs)
    write(io, " ")
    if length(terms) > 1
        write(io, "(")
        show_terms(io, terms)
        write(io, ")")
    else
        show_terms(io, terms)
    end
    return io
end

function SymbolicUtils.show_call(io::IO, ::typeof(expim), x::SymbolicUtils.BasicSymbolic; kw...)
    arg = only(SymbolicUtils.arguments(x))
    neg = leading_sign(arg) < 0
    write(io, neg ? "exp(-im*" : "exp(im*")
    body = neg ? SymbolicUtils.unwrap(expand(-arg)) : arg
    paren = SymbolicUtils.iscall(body) && SymbolicUtils.operation(body) === (+)
    paren && write(io, "(")
    print(io, body)
    paren && write(io, ")")
    write(io, ")")
    return io
end

function SymbolicUtils.show_call(io::IO, ::SumFunc, x::SymbolicUtils.BasicSymbolic; kw...)
    show_sum_prefix(io, sum_indices(x), sum_ne(x))
    write(io, " ")
    body = sum_body(x)
    paren = SymbolicUtils.iscall(body) && SymbolicUtils.operation(body) === (+)
    paren && write(io, "(")
    print(io, body)
    paren && write(io, ")")
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

show_unitary_time(::IO, ::UnitaryTransform{StaticTime}) = false
function show_unitary_time(io::IO, U::UnitaryTransform{DynamicTime})
    write(io, "time = ")
    show(io, U.time.variable)
    return true
end

function show_unitary_metadata(io::IO, U::UnitaryTransform)
    shown = show_unitary_time(io, U)
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
    limited = get(io, :limit, false) && n > SHOW_RULE_LIMIT
    if limited
        print(io, n, " rules")
    else
        for (k, generator) in enumerate(gs)
            k > 1 && write(io, ", ")
            show_unitary_rule(io, U, generator)
        end
    end
    if U isa UnitaryTransform{DynamicTime} || !iszero(U.gauge)
        write(io, "; ")
        show_unitary_metadata(io, U)
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

    limited = get(io, :limit, false) && n > SHOW_RULE_LIMIT
    if !limited
        write(io, ":")
        for generator in gs
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
