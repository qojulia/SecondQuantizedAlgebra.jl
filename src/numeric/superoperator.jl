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

# Side matrices are cached by the op vector. `Index` equality ignores the concrete site, so
# the resolved slots go into the key too; otherwise unrolled sites would share one matrix.
const SideKey = Tuple{Vector{Op}, Vector{Int}}
side_key(ops::Vector{Op}, ctx::NumericContext) = (ops, Int[resolve_slot(o, ctx.sites) for o in ops])

"""
    side_matrix(ops, ctx, cache::Dict{SideKey, M}) -> M

The matrix of the operator product `ops` in the backend of `ctx`, converted to the cache's
matrix type `M` and cached by [`side_key`](@ref).
"""
function side_matrix(ops::Vector{Op}, ctx::NumericContext, cache::Dict{SideKey, M}) where {M}
    return get!(cache, side_key(ops, ctx)) do
        op = static_operator(OpTerm[(1.0 + 0.0im, ops)], ctx, EagerOperator())
        convert(M, numeric_matrix(ctx.backend, op))
    end
end

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
[`to_numeric`](@ref) form), and `scalar_subs` resolves coefficients of unrolled indices,
such as `g(k)`. The jumps are expanded in their operator monomials `M_m`, and the rate matrix
`Γ[m, n]` of every coupled block is checked to be Hermitian and positive semidefinite and
diagonalized, `c_k = √λ_k Σ_m U[m, k] M_m`; a single channel gives `√γ J`. Eigenvalues with
`|λ| ≤ tol * max|λ|` are clipped to zero (`tol = √eps` by default), and a more negative
eigenvalue or a non-Hermitian `H` throws. A time-dependent `H` is checked at a few sample
times; with `(p, t)` coefficients it is checked symbolically, treating parameters as real.

With a `time_parameter` (QuantumToolbox only), `H` becomes a `QobjEvo` and a diagonal
channel with a time-dependent rate `γ(t)` gives the collapse operator `√γ(t) J`. Such a rate
must stay real and nonnegative: it is checked at the sample times, and evaluating a negative
or complex rate during the solve throws. A correlated channel with a time-dependent rate is
refused, because the eigenvectors of `Γ(t)` move with `t`; convert the full generator with
`to_numeric(L, ...; time_parameter)` instead.
"""
function to_numeric_master(
        L::Lindbladian{QAdd}, h::HilbertSpace, dims;
        backend::NumericBackend = default_backend(), parameter = NO_PARAMETER,
        time_parameter = NO_TIME_PARAMETER, sites::AbstractDict = EMPTY_SITES,
        scalar_subs::AbstractDict = NO_SCALAR_SUBS, tol::Real = sqrt(eps()),
    )
    basis = numeric_basis(backend, h, dims)
    site_map = Dict{Int, Vector{Int}}(sites)
    isempty(site_map) || validate_sites(backend, basis, site_map)
    rules = expand_parameter(parameter)
    tp = normalize_time_parameter(time_parameter)
    ctx = NumericContext(backend, basis, NO_SUBS, NO_SCALAR_SUBS, site_map)
    coeffs = CoefficientRules(rules, split_scalar_subs(scalar_subs)...)
    H = master_hamiltonian(L.hamiltonian, ctx, coeffs, tp)
    c_ops = master_collapse_operators(L, ctx, coeffs, tp, tol)
    return (hamiltonian = H, collapse_operators = c_ops)
end

to_numeric_master(::Lindbladian{QExpr}, h::HilbertSpace, dims; kwargs...) = formal_map_numeric_error()

# The `parameter` rules and the `scalar_subs` that resolve coefficients after unrolling.
struct CoefficientRules
    parameter::Tuple{OpRules, ScalarRules}
    sub_re::ScalarRules
    sub_im::ScalarRules
    has_imag::Bool
end

resolve_coefficient(c::CNum, sub_coef::Dict{Index, Index}, r::CoefficientRules) =
    apply_scalar_subs(change_index(c, sub_coef), r.sub_re, r.sub_im, r.has_imag)

function unrolled_terms(q::QAdd, ctx::NumericContext, r::CoefficientRules)
    terms = CoeffTerm[]
    for (t, c) in q.arguments
        dep = bound_indices(t, c, q.indices)
        for_each_assignment(dep, t.ne, ctx.sites) do sub_op, sub_coef
            nc = resolve_coefficient(c, sub_coef, r)
            push!(terms, (nc, Op[change_index(o, sub_op) for o in t.ops]))
        end
    end
    return terms
end

@noinline non_hermitian_error() = throw(
    ArgumentError("`to_numeric_master` requires a Hermitian Hamiltonian"),
)

function master_hamiltonian(H::QAdd, ctx::NumericContext, r::CoefficientRules, tp)
    H_ = apply_parameter(H, r.parameter)
    terms = unrolled_terms(H_, ctx, r)
    if !isempty(tp)
        check_td_hermitian(H_, terms, ctx, tp)
        return to_numeric_td(terms, ctx, tp)
    end
    op = static_operator(OpTerm[(to_complex(c), ops) for (c, ops) in terms], ctx, EagerOperator())
    hermitian_matrix(numeric_matrix(ctx.backend, op)) || non_hermitian_error()
    return op
end

# Sample times for the checks of a time-dependent Hamiltonian and rate: generic points, so a
# coefficient pair that agrees only at isolated instants is still caught.
const HERMITIAN_SAMPLE_TIMES = (0.0, 0.5772156649, 1.4142135624, 3.1415926536, 7.3890560989)

coefficient_at(c::ComplexF64, ::Float64) = c
coefficient_at(f, t::Float64) = ComplexF64(f(t))

# Hermitian up to `√eps` relative to `scale`.
hermitian_matrix(m, scale::Real) = maximum(abs, m - m'; init = 0.0) <= sqrt(eps()) * scale
hermitian_matrix(m) = hermitian_matrix(m, max(1.0, maximum(abs, m; init = 0.0)))

# Time values are numeric, so `H(t)` is assembled and compared with its adjoint at sample
# times. `(p, t)` coefficients cannot be evaluated before the solver runs; they fall back to
# the symbolic check, which treats parameters as real like `adjoint` does.
function check_td_hermitian(H::QAdd, terms::Vector{CoeffTerm}, ctx::NumericContext, tp)
    basevars, valuefuncs, p_aware = time_basis(tp)
    if p_aware
        iszero(simplify(H - adjoint(H))) || non_hermitian_error()
        return nothing
    end
    coeffs = [td_coeff(to_num(c), basevars, valuefuncs, false) for (c, _) in terms]
    mats = [
        numeric_matrix(ctx.backend, static_operator(OpTerm[(1.0 + 0.0im, ops)], ctx, EagerOperator()))
            for (_, ops) in terms
    ]
    isempty(mats) && return nothing
    for t in HERMITIAN_SAMPLE_TIMES
        m = sum(coefficient_at(coeffs[k], t) * mats[k] for k in eachindex(mats))
        hermitian_matrix(m) || non_hermitian_error()
    end
    return nothing
end

# Collapse operators are factored in a basis of operator monomials on concrete subsystems, so
# linearly dependent jumps such as `J` and `2J` share one rate-matrix row. A monomial's ops
# are stably sorted by subsystem: ops on different subsystems commute, and the order on one
# subsystem is physical.
const MonomialKey = Vector{NTuple{8, Int}}

monomial_op_key(o::Op, slot::Int) = (
    Int(o.kind), Int(o.name_id), Int(o.space_index), slot, Int(o.l1), Int(o.l2), Int(o.g),
    Int(o.nlev),
)

struct MonomialBasis
    index::Dict{MonomialKey, Int}
    ops::Vector{Vector{Op}}
end
MonomialBasis() = MonomialBasis(Dict{MonomialKey, Int}(), Vector{Op}[])

function monomial!(basis::MonomialBasis, ops::Vector{Op}, ctx::NumericContext)
    slots = Int[resolve_slot(o, ctx.sites) for o in ops]
    perm = sortperm(slots)
    key = MonomialKey([monomial_op_key(ops[p], slots[p]) for p in perm])
    return get!(basis.index, key) do
        push!(basis.ops, ops[perm])
        length(basis.ops)
    end
end

# A jump as its coefficients over the monomial basis.
const JumpVector = Dict{Int, ComplexF64}

function unrolled_jump(
        J::QAdd, sub_op::Dict{Index, Index}, sub_coef::Dict{Index, Index}, r::CoefficientRules,
        basis::MonomialBasis, ctx::NumericContext,
    )
    J_ = apply_parameter(J, r.parameter)
    u = JumpVector()
    for (t, c) in J_.arguments
        m = monomial!(basis, Op[change_index(o, sub_op) for o in t.ops], ctx)
        u[m] = get(u, m, 0.0im) + to_complex(resolve_coefficient(c, sub_coef, r))
    end
    return u
end

jump_terms(u::JumpVector, basis::MonomialBasis) =
    OpTerm[(c, basis.ops[m]) for (m, c) in u if !iszero(c)]

@noinline not_psd_error(λ) = throw(
    ArgumentError(
        "`to_numeric_master`: the rate matrix is not positive semidefinite " *
            "(eigenvalue $λ); a solver needs a GKSL generator",
    ),
)

@noinline symbolic_rate_error() = throw(
    ArgumentError(
        "`to_numeric_master`: a rate is symbolic; supply its value via `parameter` " *
            "or `time_parameter`",
    ),
)

@noinline non_hermitian_rates_error() = throw(
    ArgumentError("`to_numeric_master`: the rate matrix is not Hermitian"),
)

@noinline td_block_error() = throw(
    ArgumentError(
        "`to_numeric_master`: correlated channels with time-dependent rates are not " *
            "supported, because the eigenvectors of Γ(t) move with t; convert the generator " *
            "with `to_numeric(L, ...; time_parameter)` and pass the superoperator to `mesolve`",
    ),
)

# The generator is a sum over channels, so a diagonal channel with a time-dependent rate gives
# its own `√γ(t) J`, and the static channels are factored together from one rate matrix
# `Γ[m, n] = Σ γ u[m] conj(v[n])` over the monomials.
function master_collapse_operators(
        L::Lindbladian{QAdd}, ctx::NumericContext, r::CoefficientRules, tp, tol::Real,
    )
    basis = MonomialBasis()
    Γ = Dict{Tuple{Int, Int}, ComplexF64}()
    td = Tuple{JumpVector, Complex{Num}}[]
    for (ch, γ) in L.channels
        for_each_assignment(ch.scope, ch.ne, ctx.sites) do sub_op, sub_coef
            rate = apply_parameter(to_num(change_index(γ, sub_coef)), r.parameter)
            rate = to_num(apply_scalar_subs(to_cnum(rate), r.sub_re, r.sub_im, r.has_imag))
            u = unrolled_jump(ch.jump, sub_op, sub_coef, r, basis, ctx)
            if !coeff_is_const(rate)
                isempty(tp) && symbolic_rate_error()
                is_diagonal(ch) || td_block_error()
                push!(td, (u, rate))
                return nothing
            end
            v = is_diagonal(ch) ? u : unrolled_jump(ch.partner, sub_op, sub_coef, r, basis, ctx)
            γn = to_complex(rate)
            for (m, cu) in u, (n, cv) in v
                Γ[(m, n)] = get(Γ, (m, n), 0.0im) + γn * cu * conj(cv)
            end
            return nothing
        end
    end
    static_terms = Vector{OpTerm}[]
    for block in rate_blocks(length(basis.ops), Γ)
        append_block_terms!(static_terms, block, Γ, basis, tol)
    end
    static_ops = [static_operator(terms, ctx, EagerOperator()) for terms in static_terms]
    td_ops = [td_collapse(jump_terms(u, basis), rate, ctx, tp) for (u, rate) in td]
    return vcat(static_ops, td_ops)
end

function block_root(parent::Vector{Int}, x::Int)
    while parent[x] != x
        x = parent[x]
    end
    return x
end

# Groups of monomials coupled by nonzero rates; each is factored on its own.
function rate_blocks(n::Int, Γ::Dict{Tuple{Int, Int}, ComplexF64})
    parent = collect(1:n)
    for ((a, b), γ) in Γ
        iszero(γ) && continue
        ra, rb = block_root(parent, a), block_root(parent, b)
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

# Diagonalize the rate matrix of one block, `c_k = √λ_k Σ_m U[m, k] M_m`.
function append_block_terms!(
        out::Vector{Vector{OpTerm}}, block::Vector{Int}, Γ::Dict{Tuple{Int, Int}, ComplexF64},
        basis::MonomialBasis, tol::Real,
    )
    G = ComplexF64[get(Γ, (a, b), 0.0im) for a in block, b in block]
    scale_ = maximum(abs, G; init = 0.0)
    iszero(scale_) && return out
    hermitian_matrix(G, scale_) || non_hermitian_rates_error()
    F = LinearAlgebra.eigen(LinearAlgebra.Hermitian((G + G') / 2))
    cutoff = tol * maximum(abs, F.values)
    for k in eachindex(F.values)
        λ = F.values[k]
        λ < -cutoff && not_psd_error(λ)
        λ <= cutoff && continue
        terms = OpTerm[]
        for (row, m) in enumerate(block)
            w = sqrt(λ) * F.vectors[row, k]
            iszero(w) || push!(terms, (w, basis.ops[m]))
        end
        push!(out, terms)
    end
    return out
end

@noinline td_rate_error(t, r) = throw(
    ArgumentError(
        "`to_numeric_master`: the time-dependent rate is $r at t = $t; a collapse operator " *
            "√γ(t) J needs a real, nonnegative rate",
    ),
)

# `√γ(t)` of a rate that must stay real and nonnegative; roundoff below zero is clipped.
function collapse_amplitude(r::Number, t)
    tol = sqrt(eps()) * max(1.0, abs(r))
    (abs(imag(r)) <= tol && real(r) >= -tol) || td_rate_error(t, r)
    return sqrt(max(real(r), 0.0))
end

td_collapse(::Vector{OpTerm}, ::Complex{Num}, ::NumericContext, ::Tuple{}) = symbolic_rate_error()

function td_collapse(terms::Vector{OpTerm}, rate::Complex{Num}, ctx::NumericContext, tp::AbstractDict)
    op = static_operator(terms, ctx, EagerOperator())
    basevars, valuefuncs, p_aware = time_basis(tp)
    check_time_variables(rate, basevars)
    f = compile_coeff(rate, basevars...)
    if p_aware
        g = let f = f, valuefuncs = valuefuncs
            (p, t) -> collapse_amplitude(f(map(v -> v(p, t), valuefuncs)...), t)
        end
        return numeric_scaled_td(ctx.backend, op, g)
    end
    g = let f = f, valuefuncs = valuefuncs
        t -> collapse_amplitude(f(map(v -> v(t), valuefuncs)...), t)
    end
    foreach(g, HERMITIAN_SAMPLE_TIMES)
    return numeric_scaled_td(ctx.backend, op, g)
end
