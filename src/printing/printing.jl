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
Base.show(io::IO, h::PauliSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::SpinSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))
Base.show(io::IO, h::PhaseSpace) = (write(io, "ℋ("); print(io, h.name); write(io, ")"))

Base.show(io::IO, idx::Index) = print(io, idx.name)

function _show_index_suffix(io::IO, idx::Index)
    if has_index(idx)
        write(io, "_")
        print(io, idx.name)
    end
    return
end

function Base.show(io::IO, x::Destroy)
    print(io, x.name)
    return _show_index_suffix(io, x.index)
end
function Base.show(io::IO, x::Create)
    print(io, x.name)
    _show_index_suffix(io, x.index)
    return write(io, "'")
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
function Base.show(io::IO, x::Transition)
    print(io, x.name)
    _show_index_suffix(io, x.index)
    _write_subscript(io, x.i)
    return _write_subscript(io, x.j)
end

const _xyz_sym = (:x, :y, :z)
function Base.show(io::IO, x::Pauli)
    print(io, x.name)
    _show_index_suffix(io, x.index)
    return print(io, _xyz_sym[x.axis])
end
function Base.show(io::IO, x::Spin)
    print(io, x.name)
    _show_index_suffix(io, x.index)
    return print(io, _xyz_sym[x.axis])
end

function Base.show(io::IO, x::Position)
    print(io, x.name)
    return _show_index_suffix(io, x.index)
end
function Base.show(io::IO, x::Momentum)
    print(io, x.name)
    return _show_index_suffix(io, x.index)
end

function _show_prefactor(io::IO, c::CNum)
    return if iszero(imag(c))
        print(io, real(c))
    elseif iszero(real(c))
        i = imag(c)
        if isone(i)
            write(io, "im")
        elseif isequal(i, Num(-1))
            write(io, "-im")
        else
            print(io, i)
            write(io, "im")
        end
    else
        write(io, "(")
        print(io, real(c))
        write(io, " + ")
        print(io, imag(c))
        write(io, "im)")
    end
end
_show_prefactor(io::IO, c::Number) = print(io, c)

function _is_unit(c::CNum)
    return iszero(imag(c)) && isone(real(c))
end
_is_unit(c::Number) = isone(c)

function _is_neg_unit(c::CNum)
    return iszero(imag(c)) && isequal(real(c), Num(-1))
end
_is_neg_unit(c::Number) = isequal(c, -1)

function _is_real_negative(c::CNum)
    if _iszero_num(c.im)
        r = Symbolics.value(SymbolicUtils.unwrap(real(c)))
        return r isa Real && r < 0
    end
    return false
end
_is_real_negative(::Number) = false

# Check if a symbolic prefactor needs parentheses when followed by operators.
function _needs_pf_parens(c::CNum)
    iszero(imag(c)) || return false
    r = Symbolics.value(SymbolicUtils.unwrap(real(c)))
    r isa Number && return false
    SymbolicUtils.iscall(r) || return false
    op = SymbolicUtils.operation(r)
    return op === (/) || op === (+)
end

function _show_term(io::IO, c::CNum, ops::Vector{QSym})
    if isempty(ops)
        _show_prefactor(io, c)
        return
    end
    if _is_neg_unit(c)
        write(io, "-")
    elseif !_is_unit(c)
        if _needs_pf_parens(c)
            write(io, "(")
            _show_prefactor(io, c)
            write(io, ")")
        else
            _show_prefactor(io, c)
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
        print(io, idx.name)
        write(io, "=1:")
        print(io, idx.range)
        write(io, ")")
    end
    if !isempty(ne_pairs)
        write(io, "(")
        for (k, (a, b)) in enumerate(ne_pairs)
            k > 1 && write(io, ",")
            print(io, a.name)
            write(io, "≠")
            print(io, b.name)
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

function SymbolicUtils.show_metadata(
        io::IO, x::SymbolicUtils.BasicSymbolic, ::Type{SumIndices}, val::Vector{Index},
    )
    isempty(val) && return false
    _show_sum_prefix(io, val, _restore_sum_metadata_ne(x))
    write(io, " ")
    SymbolicUtils.show_plain(io, x)
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
