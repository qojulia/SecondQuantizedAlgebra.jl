# Numeric conversion of superoperators and Lindbladians, and collapse-operator solver inputs.

const SandwichTerm = Tuple{ComplexF64, Vector{Op}, Vector{Op}}

"""
    numeric_superoperator(be, ctx::NumericContext, terms::Vector{Tuple{ComplexF64, Vector{Op}, Vector{Op}}})

Build the backend's eager superoperator `Σ c * sprepost(L, R)` from sandwich terms. Each side
is converted with the ordinary operator path of `ctx`.
"""
function numeric_superoperator end

"""
    numeric_scaled_td(be, op, f)

Wrap the static backend operator `op` with the time-dependent scalar coefficient `f`, a
`t -> value` or `(p, t) -> value` closure. Only backends with time-dependent collapse
operators implement it.
"""
function numeric_scaled_td end

numeric_scaled_td(::NumericBackend, op, f) = throw(
    ArgumentError(
        "time-dependent collapse operators need the QuantumToolbox backend; " *
            "QuantumOpticsBase supports only static solver inputs",
    ),
)

"""
    numeric_matrix(be, op) -> AbstractMatrix

The matrix of an eager backend operator, used for the Hermiticity check of
[`to_numeric_master`](@ref).
"""
function numeric_matrix end

const NumericMap = Union{SuperOperator, Lindbladian{QAdd}}
const FormalMap = Union{SuperExpr, Lindbladian{QExpr}}

@noinline formal_map_numeric_error() = throw(
    ArgumentError(
        "direct numeric conversion of maps with formal operator functions is not " *
            "implemented; call `taylor(X, 0:n)` first",
    ),
)

@noinline indexed_map_error() = throw(
    ArgumentError(
        "the map carries bound summation indices; convert it with the indexed form " *
            "`to_numeric(S, basis_or_dims, sites)`",
    ),
)

@noinline map_average_error() = throw(
    ArgumentError(
        "a superoperator has no expectation value; average the Heisenberg drift " *
            "`transpose(L)(A)` instead",
    ),
)

sandwich_terms(S::SuperOperator) =
    SandwichTerm[(to_complex(c), t.left, t.right) for (t, c) in S.arguments]

apply_parameter(S::SuperOperator, (op_rules, scalar_rules)::Tuple{OpRules, ScalarRules}) =
    isempty(op_rules) && isempty(scalar_rules) ? S : substitute_split(S, op_rules, scalar_rules)

materialize_map(m, ::EagerOperator) = m
materialize_map(m, ::typeof(identity)) = m
materialize_map(m, op_type) = op_type(m)

function to_numeric_translated(S::SuperOperator, ctx::NumericContext, parameter, time_parameter, op_type)
    isempty(S.indices) || indexed_map_error()
    S_ = apply_parameter(S, parameter)
    if isempty(time_parameter)
        return materialize_map(numeric_superoperator(ctx.backend, ctx, sandwich_terms(S_)), op_type)
    end
    check_td_op_type(op_type)
    return to_numeric_td(S_, ctx, time_parameter)
end

to_numeric_translated(L::Lindbladian{QAdd}, ctx::NumericContext, parameter, time_parameter, op_type) =
    to_numeric_translated(SuperOperator(L), ctx, parameter, time_parameter, op_type)
to_numeric_translated(::FormalMap, ::NumericContext, parameter, time_parameter, op_type) =
    formal_map_numeric_error()

to_numeric_td(::SuperOperator, ::NumericContext, time_parameter) = throw(
    ArgumentError(
        "time-dependent superoperators need the QuantumToolbox backend; " *
            "QuantumOpticsBase superoperators are static",
    ),
)

function compile_td_coeffs(S::SuperOperator, time_parameter)
    basevars, valuefuncs, p_aware = time_basis(time_parameter)
    exprs = Complex{Num}[to_num(c) for (_, c) in S.arguments]
    isempty(exprs) && push!(exprs, Complex{Num}(0))
    for e in exprs
        check_time_variables(e, basevars)
    end
    return compile_coeffs(exprs, basevars, valuefuncs, p_aware)
end

function to_numeric(
        S::Union{NumericMap, FormalMap}, b::Basis;
        parameter = NO_PARAMETER, time_parameter = NO_TIME_PARAMETER,
        operators = NO_SUBS, adjoint_ops = true, op_type = EagerOperator(),
    )
    return to_numeric_kw(QuantumOpticsBackend(), S, b; parameter, time_parameter, operators, adjoint_ops, op_type)
end

to_numeric(S::NumericMap, b::Basis, d::AbstractDict{<:QSym}) =
    to_numeric_translated(S, NumericContext(QuantumOpticsBackend(), b, d), (OpRules(), ScalarRules()), NO_TIME_PARAMETER, EagerOperator())
to_numeric(::FormalMap, ::Basis, ::AbstractDict{<:QSym}) = formal_map_numeric_error()
to_numeric(S::Union{NumericMap, FormalMap}, state::QuantumState, d::AbstractDict{<:QSym}) =
    to_numeric(S, numeric_basis(state), d)

## Indexed conversion

function for_each_assignment(f::F, dep::Vector{Index}, ne::Vector{NonEqualPair}, sites::Dict{Int, Vector{Int}}) where {F}
    if isempty(dep)
        f(Dict{Index, Index}(), Dict{Index, Index}())
        return nothing
    end
    lens = Int[length(get(sites, Int(idx.space_index), Int[])) for idx in dep]
    any(iszero, lens) && throw(
        ArgumentError("no sites supplied for the space of a bound summation index"),
    )
    sub_op = Dict{Index, Index}()
    sub_coef = Dict{Index, Index}()
    for combo in 1:prod(lens)
        empty!(sub_op)
        empty!(sub_coef)
        rem = combo - 1
        for k in eachindex(dep)
            kpos = (rem % lens[k]) + 1
            rem ÷= lens[k]
            idx = dep[k]
            sub_op[idx] = Index(idx.name_id, idx.range_id, idx.space_index, Int32(kpos))
            sub_coef[idx] = Index(Int32(0), idx.range_id, idx.space_index, Int32(kpos))
        end
        violates_ne(ne, sub_op) && continue
        f(sub_op, sub_coef)
    end
    return nothing
end

function to_numeric_indexed(
        be::NumericBackend, basis, S::SuperOperator,
        sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym}, scalar_subs::AbstractDict,
    )
    validate_sites(be, basis, sites)
    site_map = Dict{Int, Vector{Int}}(sites)
    ctx = NumericContext(be, basis, d, scalar_subs, site_map)
    sub_re, sub_im, has_imag = split_scalar_subs(scalar_subs)
    terms = SandwichTerm[]
    for (t, c) in S.arguments
        dep = bound_indices(t, c, S.indices)
        for_each_assignment(dep, t.ne, site_map) do sub_op, sub_coef
            nc = apply_scalar_subs(change_index(c, sub_coef), sub_re, sub_im, has_imag)
            push!(
                terms, (
                    to_complex(nc), Op[change_index(o, sub_op) for o in t.left],
                    Op[change_index(o, sub_op) for o in t.right],
                ),
            )
        end
    end
    return numeric_superoperator(be, ctx, terms)
end

to_numeric_indexed(
    be::NumericBackend, basis, L::Lindbladian{QAdd}, sites::AbstractDict{Int, Vector{Int}},
    d::AbstractDict{<:QSym}, scalar_subs::AbstractDict,
) = to_numeric_indexed(be, basis, SuperOperator(L), sites, d, scalar_subs)

function to_numeric(
        S::NumericMap, b::Basis, sites::AbstractDict{Int, Vector{Int}},
        d::AbstractDict{<:QSym} = NO_SUBS, scalar_subs::AbstractDict = NO_SCALAR_SUBS,
    )
    return to_numeric_indexed(QuantumOpticsBackend(), b, S, sites, d, scalar_subs)
end

## Refusals

numeric_average_impl(::Union{NumericMap, FormalMap}, ::Any, ::AbstractDict{<:QSym}) = map_average_error()
numeric_average_vec(::Union{NumericMap, FormalMap}, ::Any, ::AbstractDict{<:QSym}) = map_average_error()
expect(::Union{NumericMap, FormalMap}, state, d::AbstractDict{<:QSym} = NO_SUBS) = map_average_error()

## Solver inputs

"""
    to_numeric_master(L::Lindbladian, h, dims; backend, parameter, time_parameter, sites, tol)
        -> (hamiltonian = H, collapse_operators = c_ops)

Convert a Lindbladian to the inputs of a master-equation solver, so that
`H, c_ops = to_numeric_master(L, h, dims)` feeds QuantumOptics `timeevolution.master` and
QuantumToolbox `mesolve` without rate arguments.

Channels are unrolled over their scopes (bound indices need `sites`, as in the indexed
[`to_numeric`](@ref) form). The rate matrix of every correlated block of jump operators is
checked to be Hermitian and positive semidefinite and diagonalized,
`c_k = √λ_k Σ_a U[a, k] J_a`; a diagonal channel gives `√γ J`. Eigenvalues with
`λ ≥ -tol` are clipped to zero, where `tol = √eps * max|λ|` by default, and a more negative
eigenvalue or a non-Hermitian static `H` throws. `scalar_subs` resolves coefficients of
unrolled indices, such as `g(k)`, as in the indexed `to_numeric`.

With a `time_parameter` (QuantumToolbox only), `H` becomes a `QobjEvo` and a diagonal
channel with a time-dependent rate `γ(t)` gives the collapse operator `√γ(t) J`; the user
must keep such a rate nonnegative. A time-dependent correlated block is refused, because the
eigenvectors of `Γ(t)` move with `t`; convert the full generator with
`to_numeric(L, ...; time_parameter)` instead.
"""
function to_numeric_master(
        L::Lindbladian{QAdd}, h::HilbertSpace, dims;
        backend::NumericBackend = default_backend(), parameter = NO_PARAMETER,
        time_parameter = NO_TIME_PARAMETER, sites::AbstractDict = EMPTY_SITES,
        scalar_subs::AbstractDict = NO_SCALAR_SUBS, tol::Union{Nothing, Real} = nothing,
    )
    basis = numeric_basis(backend, h, dims)
    site_map = Dict{Int, Vector{Int}}(sites)
    isempty(site_map) || validate_sites(backend, basis, site_map)
    rules = expand_parameter(parameter)
    tp = normalize_time_parameter(time_parameter)
    sub_re, sub_im, has_imag = split_scalar_subs(scalar_subs)
    ctx = NumericContext(backend, basis, NO_SUBS, NO_SCALAR_SUBS, site_map)
    H = master_hamiltonian(L.hamiltonian, ctx, rules, tp, site_map, sub_re, sub_im, has_imag)
    c_ops = master_collapse_operators(L, ctx, rules, tp, site_map, sub_re, sub_im, has_imag, tol)
    return (hamiltonian = H, collapse_operators = c_ops)
end

to_numeric_master(::Lindbladian{QExpr}, h::HilbertSpace, dims; kwargs...) = formal_map_numeric_error()

function unrolled_terms(
        q::QAdd, rules, site_map::Dict{Int, Vector{Int}},
        sub_re::ScalarRules, sub_im::ScalarRules, has_imag::Bool,
    )
    q_ = apply_parameter(q, rules)
    terms = Tuple{CNum, Vector{Op}}[]
    for (t, c) in q_.arguments
        dep = bound_indices(t, c, q_.indices)
        for_each_assignment(dep, t.ne, site_map) do sub_op, sub_coef
            nc = apply_scalar_subs(change_index(c, sub_coef), sub_re, sub_im, has_imag)
            push!(terms, (nc, Op[change_index(o, sub_op) for o in t.ops]))
        end
    end
    return terms
end

@noinline non_hermitian_error() = throw(
    ArgumentError("`to_numeric_master` requires a Hermitian Hamiltonian"),
)

function master_hamiltonian(
        H::QAdd, ctx::NumericContext, rules, tp, site_map::Dict{Int, Vector{Int}},
        sub_re::ScalarRules, sub_im::ScalarRules, has_imag::Bool,
    )
    if !isempty(tp)
        isempty(H.indices) || indexed_map_error()
        return to_numeric_td(apply_parameter(H, rules), ctx, tp)
    end
    terms = OpTerm[(to_complex(c), ops) for (c, ops) in unrolled_terms(H, rules, site_map, sub_re, sub_im, has_imag)]
    op = static_operator(terms, ctx, EagerOperator())
    m = numeric_matrix(ctx.backend, op)
    scale_ = max(1.0, maximum(abs, m; init = 0.0))
    maximum(abs, m - m'; init = 0.0) <= sqrt(eps()) * scale_ || non_hermitian_error()
    return op
end

op_site_key(o::Op) = (
    Int(o.kind), Int(o.name_id), Int(o.space_index), Int(index_slot(o.index)),
    Int(o.l1), Int(o.l2), Int(o.g), Int(o.nlev), has_index(o.index) ? 1 : 0,
)

struct UnrolledJump
    terms::Vector{OpTerm}
    key::Vector{Tuple{Vector{NTuple{9, Int}}, ComplexF64}}
end

function unrolled_jump(J::QAdd, sub_op::Dict{Index, Index}, sub_coef::Dict{Index, Index}, rules, sub_re::ScalarRules, sub_im::ScalarRules, has_imag::Bool)
    J_ = apply_parameter(J, rules)
    terms = OpTerm[]
    for (t, c) in J_.arguments
        nc = apply_scalar_subs(change_index(c, sub_coef), sub_re, sub_im, has_imag)
        push!(terms, (to_complex(nc), Op[change_index(o, sub_op) for o in t.ops]))
    end
    key = [(map(op_site_key, ops), c) for (c, ops) in terms]
    sort!(key)
    return UnrolledJump(terms, key)
end

struct RateEntry
    a::Int
    b::Int
    rate::Complex{Num}
end

function master_collapse_operators(
        L::Lindbladian{QAdd}, ctx::NumericContext, rules, tp, site_map::Dict{Int, Vector{Int}},
        sub_re::ScalarRules, sub_im::ScalarRules, has_imag::Bool, tol,
    )
    jumps = UnrolledJump[]
    keys_ = Dict{Vector{Tuple{Vector{NTuple{9, Int}}, ComplexF64}}, Int}()
    entries = RateEntry[]
    function jump_slot(j::UnrolledJump)
        slot = get(keys_, j.key, 0)
        slot == 0 || return slot
        push!(jumps, j)
        keys_[j.key] = length(jumps)
        return length(jumps)
    end
    for (ch, γ) in L.channels
        dep = Index[idx for idx in ch.scope]
        for_each_assignment(dep, ch.ne, site_map) do sub_op, sub_coef
            rate = apply_parameter(to_num(change_index(γ, sub_coef)), rules)
            rate = to_num(apply_scalar_subs(to_cnum(rate), sub_re, sub_im, has_imag))
            a = jump_slot(unrolled_jump(ch.jump, sub_op, sub_coef, rules, sub_re, sub_im, has_imag))
            b = is_diagonal(ch) ? a :
                jump_slot(unrolled_jump(ch.partner, sub_op, sub_coef, rules, sub_re, sub_im, has_imag))
            push!(entries, RateEntry(a, b, rate))
        end
    end
    return collapse_from_rates(jumps, entries, ctx, tp, tol)
end

function block_root(parent::Vector{Int}, x::Int)
    while parent[x] != x
        x = parent[x]
    end
    return x
end

function rate_blocks(n::Int, entries::Vector{RateEntry})
    parent = collect(1:n)
    for e in entries
        ra, rb = block_root(parent, e.a), block_root(parent, e.b)
        ra == rb || (parent[ra] = rb)
    end
    groups = Dict{Int, Vector{Int}}()
    for k in 1:n
        push!(get!(groups, block_root(parent, k), Int[]), k)
    end
    blocks = collect(values(groups))
    sort!(blocks; by = first)
    return blocks
end

@noinline not_psd_error(λ) = throw(
    ArgumentError(
        "`to_numeric_master`: the rate matrix is not positive semidefinite " *
            "(eigenvalue $λ); a solver needs a GKSL generator",
    ),
)

@noinline td_block_error() = throw(
    ArgumentError(
        "`to_numeric_master`: correlated channels with time-dependent rates are not " *
            "supported, because the eigenvectors of Γ(t) move with t; convert the generator " *
            "with `to_numeric(L, ...; time_parameter)` and pass the superoperator to `mesolve`",
    ),
)

rate_is_static(r::Complex{Num}) = coeff_is_const(r)

function collapse_from_rates(jumps::Vector{UnrolledJump}, entries::Vector{RateEntry}, ctx::NumericContext, tp, tol)
    c_ops = Any[]
    for block in rate_blocks(length(jumps), entries)
        pos = Dict(a => k for (k, a) in enumerate(block))
        block_entries = [e for e in entries if haskey(pos, e.a)]
        static = all(e -> rate_is_static(e.rate), block_entries)
        if !static
            isempty(tp) && throw(
                ArgumentError(
                    "`to_numeric_master`: a rate is symbolic; supply its value via `parameter` " *
                        "or `time_parameter`",
                ),
            )
            length(block) == 1 || td_block_error()
            rate = sum(e.rate for e in block_entries)
            push!(c_ops, td_collapse(jumps[only(block)], rate, ctx, tp))
            continue
        end
        Γ = zeros(ComplexF64, length(block), length(block))
        for e in block_entries
            Γ[pos[e.a], pos[e.b]] += to_complex(e.rate)
        end
        scale_ = maximum(abs, Γ; init = 0.0)
        iszero(scale_) && continue
        maximum(abs, Γ - Γ'; init = 0.0) <= sqrt(eps()) * scale_ || throw(
            ArgumentError("`to_numeric_master`: the rate matrix is not Hermitian"),
        )
        F = LinearAlgebra.eigen(LinearAlgebra.Hermitian((Γ + Γ') / 2))
        cutoff = tol === nothing ? sqrt(eps()) * maximum(abs, F.values) : Float64(tol)
        for k in eachindex(F.values)
            λ = F.values[k]
            λ < -cutoff && not_psd_error(λ)
            λ <= cutoff && continue
            terms = OpTerm[]
            for (row, a) in enumerate(block)
                w = sqrt(λ) * F.vectors[row, k]
                iszero(w) && continue
                for (c, ops) in jumps[a].terms
                    push!(terms, (w * c, ops))
                end
            end
            push!(c_ops, static_operator(terms, ctx, EagerOperator()))
        end
    end
    return [c for c in c_ops]
end

function td_collapse(j::UnrolledJump, rate::Complex{Num}, ctx::NumericContext, tp)
    op = static_operator(j.terms, ctx, EagerOperator())
    basevars, valuefuncs, p_aware = time_basis(tp)
    check_time_variables(rate, basevars)
    f = compile_coeff(rate, basevars...)
    g = p_aware ?
        let f = f, valuefuncs = valuefuncs
            (p, t) -> sqrt(real(f(map(v -> v(p, t), valuefuncs)...)))
    end :
        let f = f, valuefuncs = valuefuncs
            t -> sqrt(real(f(map(v -> v(t), valuefuncs)...)))
    end
    return numeric_scaled_td(ctx.backend, op, g)
end
