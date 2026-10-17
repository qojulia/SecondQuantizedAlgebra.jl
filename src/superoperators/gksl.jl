"""
    GKSLFormError(message)

The error thrown when a map has no GKSL form.

The map either does not preserve the trace, has a non-Hermitian Kossakowski matrix, has a
one-sided part that is not a commutator, or has a Hamiltonian that is not Hermitian modulo the
identity. The trace is checked first, so a trace defect is reported as
such. A Hamiltonian such as `x * p` is accepted, since it differs from the Hermitian
`(x * p + p * x) / 2` by a scalar.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: spre

julia> @qnumbers a::Destroy(FockSpace(:c));

julia> kossakowski(spre(a), JumpBasis(a))
ERROR: the map does not preserve the trace
[...]
```

See also [`kossakowski`](@ref), [`Lindbladian`](@ref).
"""
struct GKSLFormError <: Exception
    message::String
end

Base.showerror(io::IO, e::GKSLFormError) = print(io, e.message)

const CoeffMatrix = Matrix{CNum}

coeff_zeros(rows::Int, columns::Int) = fill(CNUM_ZERO, rows, columns)

function coeff_identity(n::Int)
    result = coeff_zeros(n, n)
    for k in 1:n
        result[k, k] = CNUM_ONE
    end
    return result
end

function canonical_operator(operator::Union{QSym, QAdd})
    result = simplify(expand_completeness(as_side(operator)))
    isempty(result.indices) || throw(
        ArgumentError("GKSL coordinates do not support operators with bound symbolic sums")
    )
    return result
end

function scalar_part(operator::QAdd)
    scalar = CNUM_ZERO
    for (term, c) in operator.arguments
        isempty(term.ops) || continue
        scalar = simplify_prefactor(scalar + c)
    end
    return scalar
end

function identity_free(canonical::QAdd)
    scalar = scalar_part(canonical)
    return iszero(scalar) ? canonical : simplify(canonical - scalar * one(canonical))
end

function jump_direction(operator::Union{QSym, QAdd})
    projected = identity_free(canonical_operator(operator))
    iszero(projected) && throw(
        ArgumentError("a jump-basis operator cannot be proportional to the identity")
    )
    return projected
end

unit_monomial(term::QTerm) = single_qadd(CNUM_ONE, term.ops, term.ne)

function pivot_row!(work::CoeffMatrix, first_row::Int, column::Int)
    for row in first_row:size(work, 1)
        value = simplify_prefactor(work[row, column])
        work[row, column] = value
        iszero(value) || return row
    end
    return 0
end

function swap_rows!(matrix::CoeffMatrix, first_row::Int, second_row::Int)
    first_row == second_row && return matrix
    for column in axes(matrix, 2)
        matrix[first_row, column], matrix[second_row, column] =
            matrix[second_row, column], matrix[first_row, column]
    end
    return matrix
end

function coeff_inverse(matrix::CoeffMatrix)
    n = LinearAlgebra.checksquare(matrix)
    left = copy(matrix)
    right = coeff_identity(n)
    for column in 1:n
        candidate = pivot_row!(left, column, column)
        iszero(candidate) && throw(ArgumentError("jump-basis pivot block is singular"))
        swap_rows!(left, column, candidate)
        swap_rows!(right, column, candidate)
        scale = inv(left[column, column])
        for trailing in 1:n
            left[column, trailing] = simplify_prefactor(left[column, trailing] * scale)
            right[column, trailing] = simplify_prefactor(right[column, trailing] * scale)
        end
        for row in 1:n
            row == column && continue
            factor = simplify_prefactor(left[row, column])
            iszero(factor) && continue
            for trailing in 1:n
                left[row, trailing] =
                    simplify_prefactor(left[row, trailing] - factor * left[column, trailing])
                right[row, trailing] =
                    simplify_prefactor(right[row, trailing] - factor * right[column, trailing])
            end
        end
    end
    return right
end

function coeff_product(left::CoeffMatrix, right::CoeffMatrix)
    size(left, 2) == size(right, 1) || throw(DimensionMismatch("matrix dimensions do not match"))
    result = coeff_zeros(size(left, 1), size(right, 2))
    for row in axes(result, 1), column in axes(result, 2)
        value = CNUM_ZERO
        for k in axes(left, 2)
            (iszero(left[row, k]) || iszero(right[k, column])) && continue
            value = value + left[row, k] * right[k, column]
        end
        result[row, column] = simplify_prefactor(value)
    end
    return result
end

function coeff_adjoint(matrix::CoeffMatrix)
    result = coeff_zeros(size(matrix, 2), size(matrix, 1))
    for row in axes(matrix, 1), column in axes(matrix, 2)
        result[column, row] = conj(matrix[row, column])
    end
    return result
end

"""
    JumpBasis(operators...)
    JumpBasis(operators::AbstractVector)

An ordered basis of jump operators ``(F_1, …, F_q)`` in which a map is written in GKSL form,

    L(ρ) = -i[H, ρ] + Σ_{μν} d_{μν} (F_μ ρ F_ν' - {F_ν' F_μ, ρ}/2),

with `d` the Kossakowski matrix returned by [`kossakowski`](@ref). Each operator is
canonicalized with [`expand_completeness`](@ref) and taken modulo the identity, since a scalar
shift of a jump changes only the Hamiltonian. The operators must be linearly independent
modulo the identity, and an `ArgumentError` is thrown otherwise, including for an operator
proportional to the identity; they need not be orthogonal, normalized or complete. Their order fixes
the order of `d`, so bases with the same span in a different order are different. An empty
basis is the basis of a coherent map.

A `JumpBasis` iterates and indexes like a vector of its canonical operators.

# Examples

```jldoctest
julia> @qnumbers a::Destroy(FockSpace(:c));

julia> JumpBasis(a + 2, a' * a)
JumpBasis with 2 operators:
  1: a
  2: a' * a
```

See also [`independent_basis`](@ref), [`support_basis`](@ref).
"""
struct JumpBasis
    operators::Vector{QAdd}
    pivot_terms::Vector{QTerm}
    pivot_inverse::CoeffMatrix
end

function reduce_against!(vector::QTermDict, echelon::Vector{Pair{QTerm, QTermDict}})
    for (pivot, row) in echelon
        entry = get(vector, pivot, CNUM_ZERO)
        iszero(entry) && continue
        factor = simplify_prefactor(entry * inv(row[pivot]))
        for (term, c) in row
            updated = simplify_prefactor(get(vector, term, CNUM_ZERO) - factor * c)
            iszero(updated) ? delete!(vector, term) : (vector[term] = updated)
        end
    end
    return vector
end

"""
    leading_terms(operators) -> (terms, independent)

Sparse row reduction of `operators` in order: the leading monomial of every operator that is
linearly independent of those before it, and which operators are.
"""
function leading_terms(operators::Vector{QAdd})
    echelon = Pair{QTerm, QTermDict}[]
    independent = falses(length(operators))
    for (k, operator) in enumerate(operators)
        vector = reduce_against!(copy(operator.arguments), echelon)
        isempty(vector) && continue
        push!(echelon, argmin(term_order_key, collect(keys(vector))) => vector)
        independent[k] = true
    end
    return QTerm[first(row) for row in echelon], independent
end

function basis_from_operators(projected::Vector{QAdd})
    pivots, independent = leading_terms(projected)
    all(independent) || throw(
        ArgumentError("jump-basis operators are linearly dependent modulo the identity")
    )
    return basis_from_pivots(projected, pivots)
end

function basis_from_pivots(operators::Vector{QAdd}, pivots::Vector{QTerm})
    sort!(pivots; by = term_order_key)
    block = CNum[get(F.arguments, term, CNUM_ZERO) for term in pivots, F in operators]
    return JumpBasis(operators, pivots, coeff_inverse(block))
end

JumpBasis(operators::AbstractVector) =
    basis_from_operators(QAdd[jump_direction(o) for o in operators])
JumpBasis(operators::Union{QSym, QAdd}...) = JumpBasis(collect(operators))
JumpBasis() = basis_from_operators(QAdd[])

Base.length(basis::JumpBasis) = length(basis.operators)
Base.isempty(basis::JumpBasis) = isempty(basis.operators)
Base.iterate(basis::JumpBasis, state...) = iterate(basis.operators, state...)
Base.eltype(::Type{JumpBasis}) = QAdd
Base.getindex(basis::JumpBasis, k::Int) = basis.operators[k]
Base.firstindex(basis::JumpBasis) = firstindex(basis.operators)
Base.lastindex(basis::JumpBasis) = lastindex(basis.operators)
Base.keys(basis::JumpBasis) = keys(basis.operators)
Base.:(==)(a::JumpBasis, b::JumpBasis) = a.operators == b.operators
Base.isequal(a::JumpBasis, b::JumpBasis) = isequal(a.operators, b.operators)
Base.hash(basis::JumpBasis, h::UInt) = hash(:JumpBasis, hash(basis.operators, h))

"""
    independent_basis(candidates::AbstractVector) -> JumpBasis

The [`JumpBasis`](@ref) of those `candidates` that are linearly independent modulo the
identity of the candidates before them, kept in order. A candidate proportional to the
identity is skipped.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: independent_basis

julia> @qnumbers a::Destroy(FockSpace(:c));

julia> independent_basis([a, 2a + 1, one(a), a' * a])
JumpBasis with 2 operators:
  1: a
  2: a' * a
```

See also [`JumpBasis`](@ref), [`support_basis`](@ref).
"""
function independent_basis(candidates::AbstractVector)
    projected = QAdd[]
    for candidate in candidates
        candidate isa Union{Op, QAdd} || throw(
            ArgumentError("a jump-basis candidate must be an operator or a polynomial, got $(typeof(candidate))")
        )
        operator = identity_free(canonical_operator(candidate))
        iszero(operator) || push!(projected, operator)
    end
    pivots, independent = leading_terms(projected)
    return basis_from_pivots(projected[independent], pivots)
end

function canonical_map(S::SuperOperator)
    isempty(S.indices) || throw(
        ArgumentError("GKSL coordinates do not support maps with bound symbolic sums")
    )
    builder = SuperBuilder()
    for (term, coefficient) in S.arguments
        left = canonical_operator(left_operator(term))
        right = canonical_operator(right_operator(term))
        for (lt, lc) in left.arguments, (rt, rc) in right.arguments
            c = simplify_prefactor(coefficient * lc * rc)
            iszero(c) && continue
            accumulate!(builder, c, unit_monomial(lt), unit_monomial(rt))
        end
    end
    return finish!(builder)
end

two_sided(left::QAdd, right::QAdd) = !isone(left) && !isone(right)

function has_two_sided_terms(canonical::SuperOperator)
    for term in keys(canonical.arguments)
        two_sided(left_operator(term), right_operator(term)) && return true
    end
    return false
end

function sandwich_pivot_matrix(canonical::SuperOperator, basis::JumpBasis)
    q = length(basis)
    result = coeff_zeros(q, q)
    index = Dict{QTerm, Int}(term => k for (k, term) in enumerate(basis.pivot_terms))
    for (term, coefficient) in canonical.arguments
        left = left_operator(term)
        right = right_operator(term)
        two_sided(left, right) || continue
        row = get(index, first(first(left.arguments)), 0)
        iszero(row) && continue
        for (rt, rc) in canonical_operator(adjoint(right)).arguments
            isempty(rt.ops) && continue
            column = get(index, rt, 0)
            iszero(column) && continue
            result[row, column] = simplify_prefactor(result[row, column] + coefficient * conj(rc))
        end
    end
    return result
end

function dissipative_map(operators::Vector{QAdd}, matrix::CoeffMatrix)
    builder = SuperBuilder()
    for row in axes(matrix, 1), column in axes(matrix, 2)
        c = matrix[row, column]
        iszero(c) && continue
        accumulate!(builder, c * SuperOperator(dissipator(operators[row], operators[column])))
    end
    return canonical_map(finish!(builder))
end

function conjugate_pair(left::CNum, right::CNum)
    iszero(simplify_prefactor(left - conj(right))) && return true
    difference = to_num(left) - conj(to_num(right))
    return iszero(Symbolics.expand(real(difference))) &&
        iszero(Symbolics.expand(imag(difference)))
end

function hermitian_part(matrix::CoeffMatrix)
    result = copy(matrix)
    for row in axes(matrix, 1), column in row:size(matrix, 2)
        conjugate_pair(matrix[row, column], matrix[column, row]) && continue
        average = simplify_prefactor((matrix[row, column] + conj(matrix[column, row])) * CNUM_HALF)
        result[row, column] = average
        result[column, row] = simplify_prefactor(conj(average))
    end
    return result
end

const AtomPowers = Tuple{Vector{RawExpression}, Vector{Rational{Int}}}

function monomial_entries(c::CNum)
    entries = Tuple{AtomPowers, CoeffScalar}[]
    tail = c.tail
    if tail isa Native
        scalar = native_scalar(c)
        iszero(scalar) || push!(entries, ((EMPTY_SYMS, EMPTY_EXPS), scalar))
    elseif tail isa Poly
        for m in tail.terms
            scalar = monomial_scalar(m)
            iszero(scalar) || push!(entries, ((m.syms, m.exps), scalar))
        end
    end
    return entries
end

scalar_magnitude(z::ComplexF64) = abs(z)
scalar_magnitude(z::GaussianRational) = hypot(Float64(real(z)), Float64(imag(z)))

function has_float_constant(node)
    SymbolicUtils.unwrap_const(node) isa Union{Float64, ComplexF64} && return true
    SymbolicUtils.iscall(node) || return false
    return any(has_float_constant, SymbolicUtils.arguments(node))
end

function has_float(c::CNum)
    tail = c.tail
    tail isa RawSymbolicCoeff && return has_float_constant(tail.expr)
    return any(entry -> last(entry) isa ComplexF64, monomial_entries(c))
end
has_float(S::SuperOperator) = any(has_float, values(S.arguments))
has_float(matrix::CoeffMatrix) = any(has_float, matrix)

function magnitude_terms(terms::Vector{Monomial{E}}) where {E}
    return Monomial{E}[
        Monomial{E}(ComplexF64(scalar_magnitude(monomial_scalar(m))), m.syms, m.exps)
            for m in terms
    ]
end

function magnitude(c::CNum)
    is_native(c) && return native(ComplexF64(scalar_magnitude(native_scalar(c))))
    is_poly(c) && return poly_coeff(Poly(magnitude_terms((c.tail::Poly).terms)))
    return CNUM_ZERO
end

magnitude(matrix::CoeffMatrix) = map(magnitude, matrix)

function magnitude(operator::QAdd)
    result = QAddBuilder()
    for (term, c) in operator.arguments
        accumulate!(result, magnitude(c) * unit_monomial(term))
    end
    return build(result)
end

function negligible(residual::CNum, references::NTuple{2, CNum}, rtol::Float64)
    iszero(residual) && return true
    iszero(rtol) && return false
    entries = monomial_entries(residual)
    isempty(entries) && return false
    scale = Dict{AtomPowers, Float64}()
    for reference in references, (key, scalar) in monomial_entries(reference)
        scale[key] = max(get(scale, key, 0.0), scalar_magnitude(scalar))
    end
    for (key, scalar) in entries
        scalar isa ComplexF64 || return false
        scalar_magnitude(scalar) <= rtol * get(scale, key, 0.0) || return false
    end
    return true
end

struct RoundoffScale
    canonical::SuperOperator
    bound::SuperOperator
    rtol::Float64
    hint::String
end

const FLOAT_HINT =
    "; the map or the basis has floating-point coefficients, so if the mismatch is roundoff " *
    "from building the map, use exact rationals such as 1//10 or pass an explicit rtol"

function roundoff_scale(
        canonical::SuperOperator, basis::JumpBasis, sandwich::CoeffMatrix, rtol::Float64,
    )
    floats = has_float(canonical) || any(F -> any(has_float, values(F.arguments)), basis.operators)
    hint = floats && iszero(rtol) ? FLOAT_HINT : ""
    (floats && rtol > 0) || return RoundoffScale(canonical, zero(SuperOperator), 0.0, hint)
    inverse = magnitude(basis.pivot_inverse)
    bound = coeff_product(coeff_product(inverse, magnitude(sandwich)), permutedims(inverse))
    operators = QAdd[magnitude(F) for F in basis.operators]
    return RoundoffScale(canonical, dissipative_map(operators, bound), rtol, hint)
end

function drop_roundoff(residual::SuperOperator, scale::RoundoffScale)
    iszero(scale.rtol) && return residual
    builder = SuperBuilder()
    for (term, c) in residual.arguments
        references = (get(scale.canonical, term, CNUM_ZERO), get(scale.bound, term, CNUM_ZERO))
        negligible(c, references, scale.rtol) && continue
        accumulate!(builder, c, left_operator(term), right_operator(term))
    end
    return finish!(builder)
end

function hermitian_modulo_identity(H::QAdd, rtol::Float64)
    Hd = canonical_operator(adjoint(H))
    for (term, c) in canonical_operator(H - Hd).arguments
        isempty(term.ops) && continue
        references = (get(H.arguments, term, CNUM_ZERO), get(Hd.arguments, term, CNUM_ZERO))
        negligible(c, references, rtol) || return false
    end
    return true
end

function residual_hamiltonian(residual::SuperOperator, scale::RoundoffScale)
    builder = QAddBuilder()
    for (term, coefficient) in residual.arguments
        left = left_operator(term)
        (isone(right_operator(term)) && !isone(left)) || continue
        accumulate!(builder, (im * coefficient) * left)
    end
    H = simplify(build(builder))
    rest = drop_roundoff(canonical_map(residual - SuperOperator(hamiltonian_action(H))), scale)
    iszero(rest) || throw(
        GKSLFormError("the one-sided part of the map is not a Hamiltonian commutator" * scale.hint)
    )
    hermitian_modulo_identity(H, scale.rtol) || throw(
        GKSLFormError("the Hamiltonian is not Hermitian modulo the identity" * scale.hint)
    )
    return H
end

function dissipative_residual(
        canonical::SuperOperator, basis::JumpBasis, matrix::CoeffMatrix, scale::RoundoffScale,
    )
    return drop_roundoff(
        canonical_map(canonical - dissipative_map(basis.operators, matrix)), scale,
    )
end

function trace_preserving(canonical::SuperOperator, rtol::Float64)
    defects = QAddBuilder()
    bounds = QAddBuilder()
    for (term, c) in canonical.arguments
        product = canonical_operator(right_operator(term) * left_operator(term))
        accumulate!(defects, c * product)
        iszero(rtol) || accumulate!(bounds, magnitude(c) * magnitude(product))
    end
    defect = simplify(build(defects))
    bound = build(bounds)
    for (term, c) in defect.arguments
        negligible(c, (get(bound.arguments, term, CNUM_ZERO), CNUM_ZERO), rtol) || return false
    end
    return true
end

function extract_gksl(canonical::SuperOperator, basis::JumpBasis, rtol::Float64)
    rtol >= 0 || throw(ArgumentError("rtol must be nonnegative, got $rtol"))
    sandwich = sandwich_pivot_matrix(canonical, basis)
    raw = coeff_product(
        coeff_product(basis.pivot_inverse, sandwich), coeff_adjoint(basis.pivot_inverse),
    )
    scale = roundoff_scale(canonical, basis, sandwich, rtol)
    trace_preserving(canonical, rtol) ||
        throw(GKSLFormError("the map does not preserve the trace" * scale.hint))
    matrix = hermitian_part(raw)
    residual = dissipative_residual(canonical, basis, matrix, scale)
    if has_two_sided_terms(residual)
        outside = isequal(matrix, raw) ||
            has_two_sided_terms(dissipative_residual(canonical, basis, raw, scale))
        outside && throw(
            ArgumentError("the map has jump directions outside the jump basis" * scale.hint)
        )
        throw(GKSLFormError("the Kossakowski matrix is not Hermitian" * scale.hint))
    end
    return residual_hamiltonian(residual, scale), matrix
end

function support_basis_canonical(canonical::SuperOperator)
    found = Set{QTerm}()
    for term in keys(canonical.arguments)
        left = left_operator(term)
        right = right_operator(term)
        two_sided(left, right) || continue
        push!(found, first(first(left.arguments)))
        for rt in keys(canonical_operator(adjoint(right)).arguments)
            isempty(rt.ops) || push!(found, rt)
        end
    end
    terms = sort!(collect(found); by = term_order_key)
    return basis_from_operators(QAdd[unit_monomial(t) for t in terms])
end

"""
    support_basis(S::SuperOperator) -> JumpBasis
    support_basis(L::Lindbladian) -> JumpBasis

The [`JumpBasis`](@ref) of the operator monomials on the two-sided sandwiches of a map, after
[`expand_completeness`](@ref), in [`term_order_key`](@ref) order. Every jump direction of the
map lies in its span, so [`Lindbladian`](@ref)`(S)` reads `S` in this basis.

# Examples

```jldoctest
julia> using SecondQuantizedAlgebra: support_basis

julia> @qnumbers b::Destroy(FockSpace(:c));

julia> support_basis(dissipator(b) + dissipator(b' * b))
JumpBasis with 2 operators:
  1: b
  2: b' * b
```

See also [`JumpBasis`](@ref), [`independent_basis`](@ref).
"""
support_basis(S::SuperOperator) = support_basis_canonical(canonical_map(S))
support_basis(L::Lindbladian{QAdd}) = support_basis(SuperOperator(L))

"""
    kossakowski(S::SuperOperator, basis::JumpBasis; rtol = 0) -> Matrix{Coeff}
    kossakowski(L::Lindbladian, basis::JumpBasis; rtol = 0) -> Matrix{Coeff}

The Hermitian Kossakowski matrix `d` of a map in the ordered `basis`:

    S(ρ) = -i[H, ρ] + Σ_{μν} d_{μν} (F_μ ρ F_ν' - {F_ν' F_μ, ρ}/2).

The extraction works in the operator algebra and introduces no Hilbert-space cutoff. It
throws an `ArgumentError` when the map has a jump direction outside `basis`, and a
[`GKSLFormError`](@ref) when the map is not of GKSL form. A `Lindbladian` is read through its
`SuperOperator`, so its own channels need not lie in `basis`.

The extraction is exact. Floating-point coefficients carry the roundoff of building the map,
so a float map is usually of GKSL form only up to that roundoff and is rejected; write rates as
exact rationals such as `1//10` or as symbolic parameters instead. A positive `rtol` opts into
a heuristic: a leftover float coefficient counts as roundoff when, for every parameter
monomial, it is at most `rtol` times the size of the coefficients it was computed from. No
single `rtol` suits every basis, since roundoff grows as the basis operators approach linear
dependence, and genuine terms below `rtol` are dropped. Exact coefficients are always compared
exactly.

# Examples

```jldoctest
julia> @qnumbers a::Destroy(FockSpace(:c));

julia> using Symbolics: @variables

julia> @variables κ::Real γ::Real;

julia> kossakowski(κ * dissipator(a) + γ * dissipator(a' * a), JumpBasis(a, a + a' * a))
2×2 Matrix{SecondQuantizedAlgebra.Coeff}:
 γ + κ  -γ
 -γ     γ
```

See also [`Lindbladian`](@ref), [`support_basis`](@ref).
"""
function kossakowski(S::SuperOperator, basis::JumpBasis; rtol::Real = 0)
    _, matrix = extract_gksl(canonical_map(S), basis, Float64(rtol))
    return matrix
end

kossakowski(L::Lindbladian{QAdd}, basis::JumpBasis; rtol::Real = 0) =
    kossakowski(SuperOperator(L), basis; rtol)

function gksl_lindbladian(canonical::SuperOperator, basis::JumpBasis, rtol::Real)
    H, matrix = extract_gksl(canonical, basis, Float64(rtol))
    builder = LindbladianBuilder{QAdd}()
    accumulate!(builder, hamiltonian_action(H))
    for row in axes(matrix, 1), column in axes(matrix, 2)
        iszero(matrix[row, column]) && continue
        accumulate!(
            builder, dissipator(basis.operators[row], basis.operators[column]),
            matrix[row, column],
        )
    end
    return finish!(builder)
end

"""
    Lindbladian(S::SuperOperator[, basis::JumpBasis]; rtol = 0)

A map in GKSL form, read back as a structured generator. The result has the Hamiltonian of
`S` and one channel `d_{μν} 𝒟[F_μ, F_ν]` for every nonzero entry of the Kossakowski matrix
in `basis`, so `SuperOperator(Lindbladian(S, basis))` equals `S` up to completeness expansion
(and up to the terms a positive `rtol` drops).
Without a basis the [`support_basis`](@ref) of `S` is used. The Hamiltonian is defined modulo
the identity, which `hamiltonian_action` drops.

Throws and handles floating-point coefficients like [`kossakowski`](@ref). The rates are not
checked for positivity.

# Examples

```jldoctest
julia> @qnumbers a::Destroy(FockSpace(:c));

julia> using SecondQuantizedAlgebra: SuperOperator

julia> Lindbladian(SuperOperator(dissipator(a + 3)))
-i[3//2im * a + -3//2im * a', ·] + 𝒟[a]
```

See also [`kossakowski`](@ref), [`support_basis`](@ref).
"""
Lindbladian(S::SuperOperator, basis::JumpBasis; rtol::Real = 0) =
    gksl_lindbladian(canonical_map(S), basis, rtol)

function Lindbladian(S::SuperOperator; rtol::Real = 0)
    canonical = canonical_map(S)
    return gksl_lindbladian(canonical, support_basis_canonical(canonical), rtol)
end

function Base.show(io::IO, basis::JumpBasis)
    print(io, "JumpBasis(")
    join(io, basis.operators, ", ")
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", basis::JumpBasis)
    n = length(basis)
    print(io, "JumpBasis with ", n, n == 1 ? " operator" : " operators")
    isempty(basis) && return nothing
    print(io, ":")
    for (k, operator) in enumerate(basis.operators)
        print(io, "\n  ", k, ": ", operator)
    end
    return nothing
end
