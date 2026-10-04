"""
    QField

Abstract supertype for all second-quantized operator expressions.

Expression layers:
- [`QSym`](@ref): atomic operator leaves; [`Op`](@ref) is the sole concrete subtype.
- [`QAdd`](@ref): canonical polynomial expressions, stored as sums of [`QTerm`](@ref)
  operator products. Ordinary polynomial arithmetic stays on this eager hot path.
- [`QExpr`](@ref): cold-path formal non-polynomial expressions such as `sin(A)`, `cos(A)`,
  and `expim(A)`. Formal functions are not expanded implicitly; use [`taylor`](@ref) to
  lower them explicitly to a finite `QAdd` polynomial.

Supports arithmetic (`+`, `-`, `*`, `^`, `/`), `adjoint`, and comparison via
`==`/`isequal`. Polynomial multiplication eagerly applies the canonical operator algebra;
operations involving a formal function preserve a `QExpr` outer representation until an
explicit lowering operation is requested.
"""
abstract type QField end

"""
    QSym <: QField

Abstract type for fundamental operator leaves.

[`Op`](@ref) is the sole concrete subtype. Its compact tagged representation stores the
operator role (`Destroy`, `Create`, `Transition`, `Pauli`, `Spin`, `Position`, `Momentum`,
etc.), display name, product-space slot, and optional symbolic index without introducing a
runtime subtype hierarchy for individual operator roles.
"""
abstract type QSym <: QField end

Base.zero(::T) where {T <: QField} = zero(T)
Base.zero(::Type{<:QField}) = 0

"""
    SiteCmp

Three-way site comparison for the canonical partial-sort. `Equal` and
`Undetermined` both mean "do not reorder," but distinguishing them in the type
prevents the "I returned 0 but meant Equal" bug class.
"""
@enum SiteCmp::UInt8 Less Equal Undetermined Greater

"""
    ReduceKind

Outcome tag for `reduce_pair`:
- `NoReduction`: pair does not reduce; `op` and `factor` slots are ignored.
- `ScalarReduction`: pair contracts to the scalar `factor`; both ops disappear.
- `OpReduction`: pair reduces to `op * factor`; `op` replaces the first input.
"""
@enum ReduceKind::UInt8 NoReduction ScalarReduction OpReduction
