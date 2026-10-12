"""
    GKSLFormError(message)

Thrown when a map has no exact GKSL coordinates: its one-sided part is not a Hamiltonian
commutator, the Kossakowski matrix is not Hermitian, or the frame is structurally singular.
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

simplify_coeff(c::CNum) = simplify_prefactor(c)

function simplify_coeffs!(matrix::CoeffMatrix)
    for k in eachindex(matrix)
        matrix[k] = simplify_coeff(matrix[k])
    end
    return matrix
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
        scalar = simplify_coeff(add_cnum(scalar, c))
    end
    return scalar
end

function traceless_direction(operator::Union{QSym, QAdd})
    canonical = canonical_operator(operator)
    scalar = scalar_part(canonical)
    projected = iszero(scalar) ? canonical : simplify(canonical - scalar * one(canonical))
    iszero(projected) && throw(
        ArgumentError("a dissipative-frame direction cannot be proportional to the identity")
    )
    return projected
end

function sorted_terms(operator::QAdd)
    entries = collect(operator.arguments)
    sort!(entries; by = entry -> term_order_key(first(entry)))
    return entries
end

function unit_monomial(term::QTerm)
    return single_qadd(CNUM_ONE, term.ops, term.ne)
end

function frame_coordinates(operators::Vector{QAdd})
    monomials = QTerm[]
    for operator in operators, (term, _) in sorted_terms(operator)
        isempty(term.ops) && continue
        term in monomials || push!(monomials, term)
    end
    coordinates = coeff_zeros(length(monomials), length(operators))
    row = Dict{QTerm, Int}(term => k for (k, term) in enumerate(monomials))
    for (column, operator) in enumerate(operators), (term, c) in operator.arguments
        isempty(term.ops) && continue
        coordinates[row[term], column] = c
    end
    return monomials, coordinates
end

function pivot_row!(work::CoeffMatrix, first_row::Int, column::Int)
    for row in first_row:size(work, 1)
        value = simplify_coeff(work[row, column])
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

function eliminate_below!(work::CoeffMatrix, pivot_row::Int, column::Int)
    pivot = work[pivot_row, column]
    for row in (pivot_row + 1):size(work, 1)
        entry = simplify_coeff(work[row, column])
        iszero(entry) && continue
        factor = simplify_coeff(entry * inv(pivot))
        for trailing in column:size(work, 2)
            work[row, trailing] = simplify_coeff(
                work[row, trailing] - factor * work[pivot_row, trailing]
            )
        end
    end
    return work
end

function pivot_monomials(coordinates::CoeffMatrix)
    monomial_count, direction_count = size(coordinates)
    direction_count == 0 && return Int[]
    monomial_count < direction_count && return Int[]
    work = permutedims(coordinates)
    pivots = Int[]
    row = 1
    for column in 1:monomial_count
        candidate = pivot_row!(work, row, column)
        iszero(candidate) && continue
        swap_rows!(work, row, candidate)
        eliminate_below!(work, row, column)
        push!(pivots, column)
        row += 1
        row > direction_count && break
    end
    return pivots
end

function columns_independent(coordinates::CoeffMatrix)
    return length(pivot_monomials(coordinates)) == size(coordinates, 2)
end

function coeff_inverse(matrix::CoeffMatrix)
    n = LinearAlgebra.checksquare(matrix)
    left = copy(matrix)
    right = coeff_identity(n)
    for column in 1:n
        candidate = pivot_row!(left, column, column)
        iszero(candidate) &&
            throw(GKSLFormError("dissipative-frame coordinate pivot is structurally singular"))
        swap_rows!(left, column, candidate)
        swap_rows!(right, column, candidate)
        scale = inv(left[column, column])
        for trailing in 1:n
            left[column, trailing] = simplify_coeff(left[column, trailing] * scale)
            right[column, trailing] = simplify_coeff(right[column, trailing] * scale)
        end
        for row in 1:n
            row == column && continue
            factor = simplify_coeff(left[row, column])
            iszero(factor) && continue
            for trailing in 1:n
                left[row, trailing] =
                    simplify_coeff(left[row, trailing] - factor * left[column, trailing])
                right[row, trailing] =
                    simplify_coeff(right[row, trailing] - factor * right[column, trailing])
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
        result[row, column] = simplify_coeff(value)
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
    DissipativeFrame(directions...)
    DissipativeFrame(directions::AbstractVector)

An ordered frame of jump directions ``(F_1, …, F_q)`` for Kossakowski coordinates,

    L(ρ) = -i[H, ρ] + Σ_{μν} d_{μν} (F_μ ρ F_ν' - {F_ν' F_μ, ρ}/2).

Each direction is canonicalized with [`expand_completeness`](@ref) and taken modulo the
identity, because a scalar shift of a jump changes only the Hamiltonian. The directions must
be linearly independent modulo the identity; they need not be orthogonal or normalized. Their
order fixes the coordinates of [`kossakowski`](@ref), so frames with the same span in a
different order are different. An empty frame is the frame of a coherent map. Read the
canonical directions with [`directions`](@ref).

See also [`independent_frame`](@ref), [`support_frame`](@ref).
"""
struct DissipativeFrame
    directions::Vector{QAdd}
    monomials::Vector{QTerm}
    coordinates::CoeffMatrix
    pivots::Vector{Int}
    pivot_inverse::CoeffMatrix
end

function frame_from_directions(projected::Vector{QAdd})
    monomials, coordinates = frame_coordinates(projected)
    pivots = pivot_monomials(coordinates)
    length(pivots) == length(projected) || throw(
        ArgumentError("dissipative-frame directions are linearly dependent modulo the identity")
    )
    pivot_matrix = coordinates[pivots, :]
    return DissipativeFrame(projected, monomials, coordinates, pivots, coeff_inverse(pivot_matrix))
end

DissipativeFrame(operators::AbstractVector) =
    frame_from_directions(QAdd[traceless_direction(o) for o in operators])
DissipativeFrame(operators::Union{QSym, QAdd}...) = DissipativeFrame(collect(operators))
DissipativeFrame() = frame_from_directions(QAdd[])

"""
    directions(frame::DissipativeFrame) -> Vector{QAdd}

The canonical directions of `frame`, in order: completeness expanded and with the identity
component removed.
"""
directions(frame::DissipativeFrame) = copy(frame.directions)

Base.length(frame::DissipativeFrame) = length(frame.directions)
Base.isempty(frame::DissipativeFrame) = isempty(frame.directions)
Base.:(==)(a::DissipativeFrame, b::DissipativeFrame) = a.directions == b.directions
Base.isequal(a::DissipativeFrame, b::DissipativeFrame) = isequal(a.directions, b.directions)
Base.hash(frame::DissipativeFrame, h::UInt) = hash(:DissipativeFrame, hash(frame.directions, h))

"""
    independent_frame(candidates) -> DissipativeFrame

The frame of the `candidates` that are linearly independent modulo the identity of the
candidates before them, kept in order. A candidate proportional to the identity, or with a
zero dissipator, is skipped.
"""
function independent_frame(candidates)
    kept = QAdd[]
    for candidate in candidates
        canonical = canonical_operator(candidate)
        projected = simplify(canonical - scalar_part(canonical) * one(canonical))
        iszero(projected) && continue
        push!(kept, projected)
        columns_independent(last(frame_coordinates(kept))) || pop!(kept)
    end
    return frame_from_directions(kept)
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
            c = simplify_coeff(coefficient * lc * rc)
            iszero(c) && continue
            accumulate!(builder, c, unit_monomial(lt), unit_monomial(rt))
        end
    end
    return build(builder)
end

two_sided(left::QAdd, right::QAdd) = !isone(left) && !isone(right)

function sandwich_pivot_matrix(canonical::SuperOperator, frame::DissipativeFrame)
    q = length(frame)
    result = coeff_zeros(q, q)
    index = Dict{QTerm, Int}(frame.monomials[p] => k for (k, p) in enumerate(frame.pivots))
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
            result[row, column] = simplify_coeff(result[row, column] + coefficient * conj(rc))
        end
    end
    return result
end

function dissipative_map(frame::DissipativeFrame, matrix::CoeffMatrix)
    q = length(frame)
    size(matrix) == (q, q) || throw(DimensionMismatch("Kossakowski matrix does not match frame"))
    builder = SuperBuilder()
    for row in 1:q, column in 1:q
        c = simplify_coeff(matrix[row, column])
        iszero(c) && continue
        accumulate!(
            builder, c * SuperOperator(dissipator(frame.directions[row], frame.directions[column])),
        )
    end
    return canonical_map(build(builder))
end

const FLOAT_ROUNDOFF_TOLERANCE = 1.0e-10

const ExactConstant = Union{Int, Rational{Int}, Complex{Int}, Complex{Rational{Int}}}
const FloatConstant = Union{Float64, ComplexF64}

constant_magnitude(value::ExactConstant) = abs(ComplexF64(value)), false
constant_magnitude(value::FloatConstant) = abs(value), true

function monomial_constant(term)
    value = SymbolicUtils.unwrap_const(term)
    value isa Union{ExactConstant, FloatConstant} && return constant_magnitude(value)
    SymbolicUtils.iscall(term) || return 1.0, false
    SymbolicUtils.operation(term) === (/) &&
        return monomial_constant(first(SymbolicUtils.arguments(term)))
    SymbolicUtils.operation(term) === (*) || return 1.0, false
    magnitude, inexact = 1.0, false
    for factor in SymbolicUtils.arguments(term)
        value = SymbolicUtils.unwrap_const(factor)
        value isa Union{ExactConstant, FloatConstant} || continue
        factor_magnitude, factor_inexact = constant_magnitude(value)
        magnitude *= factor_magnitude
        inexact |= factor_inexact
    end
    return magnitude, inexact
end

function monomial_constants!(constants::Vector{Tuple{Float64, Bool}}, part::Num)
    expanded = Symbolics.unwrap(Symbolics.expand(part))
    terms = if SymbolicUtils.iscall(expanded) && SymbolicUtils.operation(expanded) === (+)
        SymbolicUtils.arguments(expanded)
    else
        [expanded]
    end
    for term in terms
        value = SymbolicUtils.unwrap_const(term)
        value isa Union{ExactConstant, FloatConstant} && iszero(value) && continue
        push!(constants, monomial_constant(term))
    end
    return constants
end

function monomial_constants(value::CNum)
    number = to_num(value)
    constants = Tuple{Float64, Bool}[]
    monomial_constants!(constants, real(number))
    monomial_constants!(constants, imag(number))
    return constants
end

function has_floating_constant(node)
    SymbolicUtils.unwrap_const(node) isa FloatConstant && return true
    SymbolicUtils.iscall(node) || return false
    return any(has_floating_constant, SymbolicUtils.arguments(node))
end

function has_floating_constant(value::CNum)
    number = to_num(value)
    return has_floating_constant(Symbolics.unwrap(real(number))) ||
        has_floating_constant(Symbolics.unwrap(imag(number)))
end

function coefficient_scale(values)
    scale = 0.0
    for value in values, (magnitude, _) in monomial_constants(value)
        scale = max(scale, magnitude)
    end
    return iszero(scale) ? 1.0 : scale
end

function floating_roundoff(value::CNum, reference, floor::Float64)
    has_floating_constant(value) || return false
    constants = monomial_constants(value)
    isempty(constants) && return false
    threshold = FLOAT_ROUNDOFF_TOLERANCE * max(coefficient_scale(reference), floor)
    for (magnitude, inexact) in constants
        (inexact && magnitude <= threshold) || return false
    end
    return true
end

function hermitize_floating_roundoff!(matrix::CoeffMatrix)
    for row in axes(matrix, 1), column in row:size(matrix, 2)
        pair = (matrix[row, column], matrix[column, row])
        difference = simplify_coeff(matrix[row, column] - conj(matrix[column, row]))
        (iszero(difference) || !floating_roundoff(difference, pair, 0.0)) && continue
        average = simplify_coeff((matrix[row, column] + conj(matrix[column, row])) / 2)
        matrix[row, column] = average
        matrix[column, row] = simplify_coeff(conj(average))
    end
    return matrix
end

function conjugate_pair(left::CNum, right::CNum)
    difference = to_num(left) - conj(to_num(right))
    return iszero(Symbolics.expand(real(difference))) &&
        iszero(Symbolics.expand(imag(difference)))
end

function coeffs_hermitian(matrix::CoeffMatrix)
    size(matrix, 1) == size(matrix, 2) || return false
    for row in axes(matrix, 1), column in row:size(matrix, 2)
        conjugate_pair(matrix[row, column], matrix[column, row]) && continue
        pair = (matrix[row, column], matrix[column, row])
        difference = simplify_coeff(matrix[row, column] - conj(matrix[column, row]))
        iszero(difference) || floating_roundoff(difference, pair, 0.0) || return false
    end
    return true
end

function drop_floating_roundoff(
        residual::SuperOperator, canonical::SuperOperator, dissipative::SuperOperator,
        matrix::CoeffMatrix,
    )
    builder = SuperBuilder()
    matrix_scale = any(has_floating_constant, matrix) ? coefficient_scale(matrix) : 0.0
    for (term, coefficient) in residual.arguments
        left = left_operator(term)
        right = right_operator(term)
        reference = (get(canonical, term, CNUM_ZERO), get(dissipative, term, CNUM_ZERO))
        floor = two_sided(left, right) ? matrix_scale : 0.0
        floating_roundoff(coefficient, reference, floor) && continue
        accumulate!(builder, coefficient, left, right)
    end
    return build(builder)
end

function has_two_sided_terms(canonical::SuperOperator)
    for term in keys(canonical.arguments)
        two_sided(left_operator(term), right_operator(term)) && return true
    end
    return false
end

function support_terms(canonical::SuperOperator)
    found = QTerm[]
    for term in keys(canonical.arguments)
        left = left_operator(term)
        right = right_operator(term)
        two_sided(left, right) || continue
        lt = first(first(left.arguments))
        lt in found || push!(found, lt)
        for (rt, _) in canonical_operator(adjoint(right)).arguments
            isempty(rt.ops) && continue
            rt in found || push!(found, rt)
        end
    end
    sort!(found; by = term_order_key)
    return found
end

function residual_hamiltonian(canonical::SuperOperator)
    H = zero_qadd()
    for (term, coefficient) in canonical.arguments
        left = left_operator(term)
        right = right_operator(term)
        if !isone(left) && isone(right)
            H = H + (im * coefficient) * left
        elseif !(isone(left) && !isone(right))
            throw(GKSLFormError("the one-sided part of the map is not a Hamiltonian commutator"))
        end
    end
    H = simplify(H)
    iszero(canonical_map(canonical - SuperOperator(hamiltonian_action(H)))) ||
        throw(GKSLFormError("the one-sided part of the map is not a Hamiltonian commutator"))
    iszero(simplify(canonical_operator(H - adjoint(H)))) ||
        throw(GKSLFormError("the extracted Hamiltonian is not Hermitian modulo the identity"))
    return H
end

function canonical_gksl_coordinates(canonical::SuperOperator, frame::DissipativeFrame)
    sandwich = sandwich_pivot_matrix(canonical, frame)
    matrix = coeff_product(
        coeff_product(frame.pivot_inverse, sandwich), coeff_adjoint(frame.pivot_inverse),
    )
    simplify_coeffs!(matrix)
    hermitize_floating_roundoff!(matrix)
    coeffs_hermitian(matrix) || throw(
        GKSLFormError("the Kossakowski matrix is not Hermitian in the supplied frame")
    )
    dissipative = dissipative_map(frame, matrix)
    residual = canonical_map(canonical - dissipative)
    residual = drop_floating_roundoff(residual, canonical, dissipative, matrix)
    has_two_sided_terms(residual) && throw(
        ArgumentError("the map has dissipative directions outside the supplied frame")
    )
    return residual_hamiltonian(residual), matrix
end

support_frame_canonical(canonical::SuperOperator) =
    frame_from_directions(QAdd[unit_monomial(t) for t in support_terms(canonical)])

"""
    support_frame(S::SuperOperator) -> DissipativeFrame

The frame of the operator monomials that appear on the two-sided sandwiches of `S`, after
completeness expansion, in [`term_order_key`](@ref) order. Every dissipative direction of `S`
lies in its span.
"""
support_frame(S::SuperOperator) = support_frame_canonical(canonical_map(S))
support_frame(L::Lindbladian{QAdd}) = support_frame(SuperOperator(L))

"""
    gksl_coordinates(S::SuperOperator, frame::DissipativeFrame) -> (H, d)

The Hamiltonian `H` and the Kossakowski matrix `d` of `S` in `frame`, as returned separately
by `hamiltonian(Lindbladian(S, frame))` and [`kossakowski`](@ref)`(S, frame)`, from one
extraction. `H` is Hermitian modulo the identity and keeps the identity component it was read
with. Throws like [`kossakowski`](@ref).
"""
function gksl_coordinates(S::SuperOperator, frame::DissipativeFrame)
    return canonical_gksl_coordinates(canonical_map(S), frame)
end
gksl_coordinates(L::Lindbladian{QAdd}, frame::DissipativeFrame) =
    gksl_coordinates(SuperOperator(L), frame)

"""
    kossakowski(S::SuperOperator, frame::DissipativeFrame) -> Matrix{CNum}
    kossakowski(L::Lindbladian, frame::DissipativeFrame) -> Matrix{CNum}

The Hermitian Kossakowski matrix `d` of a map in the ordered `frame`, exact in the operator
algebra:

    S(ρ) = -i[H, ρ] + Σ_{μν} d_{μν} (F_μ ρ F_ν' - {F_ν' F_μ, ρ}/2).

The matrix is read from the two-sided sandwiches, and the one-sided remainder is checked to be
a Hamiltonian commutator. No Hilbert-space matrix or cutoff is introduced. Throws an
`ArgumentError` when `S` has a dissipative direction outside `frame` and a
[`GKSLFormError`](@ref) when `S` is not of GKSL form. A `Lindbladian` is read through its
`SuperOperator`, so its own channels need not match `frame`.
"""
function kossakowski(S::SuperOperator, frame::DissipativeFrame)
    _, matrix = gksl_coordinates(S, frame)
    return matrix
end
kossakowski(L::Lindbladian{QAdd}, frame::DissipativeFrame) =
    kossakowski(SuperOperator(L), frame)

function gksl_lindbladian(H::QAdd, frame::DissipativeFrame, matrix::CoeffMatrix)
    builder = LindbladianBuilder{QAdd}()
    accumulate!(builder, hamiltonian_action(H))
    for row in axes(matrix, 1), column in axes(matrix, 2)
        iszero(matrix[row, column]) && continue
        accumulate!(
            builder, dissipator(frame.directions[row], frame.directions[column]),
            matrix[row, column],
        )
    end
    return build(builder)
end

"""
    Lindbladian(S::SuperOperator[, frame::DissipativeFrame])

Read a map in GKSL form back as a structured generator. The result has the Hamiltonian of
`S` and one channel `d_{μν} 𝒟[F_μ, F_ν]` for every nonzero entry of
[`kossakowski`](@ref)`(S, frame)`, so `SuperOperator(Lindbladian(S, frame))` equals `S` up to
completeness expansion. Without a frame the [`support_frame`](@ref) of `S` is used. The
Hamiltonian is defined modulo the identity, which `hamiltonian_action` drops.

Throws like [`kossakowski`](@ref). The rates are exact but not checked for positivity.
"""
function Lindbladian(S::SuperOperator, frame::DissipativeFrame)
    H, matrix = canonical_gksl_coordinates(canonical_map(S), frame)
    return gksl_lindbladian(H, frame, matrix)
end

function Lindbladian(S::SuperOperator)
    canonical = canonical_map(S)
    has_two_sided_terms(canonical) || return hamiltonian_action(residual_hamiltonian(canonical))
    frame = support_frame_canonical(canonical)
    H, matrix = canonical_gksl_coordinates(canonical, frame)
    return gksl_lindbladian(H, frame, matrix)
end

function Base.show(io::IO, frame::DissipativeFrame)
    print(io, "DissipativeFrame(")
    join(io, frame.directions, ", ")
    return print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", frame::DissipativeFrame)
    n = length(frame)
    print(io, "DissipativeFrame with ", n, n == 1 ? " direction" : " directions")
    isempty(frame) && return nothing
    print(io, ":")
    for (k, direction) in enumerate(frame.directions)
        print(io, "\n  ", k, ": ", direction)
    end
    return nothing
end
