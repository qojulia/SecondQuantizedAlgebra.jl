using Test
using SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables, Num
using SymbolicUtils: SymbolicUtils
using Latexify: latexify
import SecondQuantizedAlgebra: Coeff, _to_cnum, to_num, _is_native, _is_poly, _conj_cnum,
    _mul_cnum, _CNUM_ONE, _recognize, _is_atom, _fold_const, _display_coeff, qadjoint,
    expim, _expim, _sup_norm

# A coefficient atom is a symbolic head the coefficient layer treats as one indivisible
# factor. Getting one right means satisfying a contract spread over ~24 sites in cnum.jl,
# monomial.jl, reduce.jl, numeric/coeff.jl and both renderers; `expim` originally satisfied
# about a third of them, and each miss was silent. This testset is that contract, executable.
#
# To add an atom: append an entry to `ATOMS` and make the testset pass. See the
# "Coefficient atoms" checklist in docs/src/devdocs.md.

@variables _u _v

# `build`      : argument -> the public coefficient
# `raw`        : argument -> the bare interned `BasicSymbolic`
# `unimodular` : does `c * conj(c) == 1` hold identically?
# `bounded`    : is |atom| <= 1 for all real arguments (so `_sup_norm` can bound it)?
const ATOMS = (
    (name = "expim", build = expim, raw = _expim, unimodular = true, bounded = true),
)

@testset "coefficient atom conformance" begin
    for atom in ATOMS
        @testset "$(atom.name)" begin
            arg = _u * _v
            c = atom.build(arg)
            b = atom.raw(SymbolicUtils.unwrap(arg))

            @testset "is a coefficient, on the polynomial tier" begin
                # A `Num` is declared real, so `conj` and `Complex * Num` would silently
                # split a complex-valued atom into unrelated real/imag halves.
                @test c isa Coeff
                @test _is_poly(c)
                @test _is_atom(b)
            end

            @testset "atom identity survives a rebuild" begin
                # `type` and `shape` take part in hash-consing, so a head with no registered
                # `promote_symtype`/`promote_shape` comes back from `maketerm` as a different
                # object. Every `substitute` and `Postwalk` rebuilds through `maketerm`.
                rebuilt = SymbolicUtils.maketerm(
                    typeof(b), SymbolicUtils.operation(b), SymbolicUtils.arguments(b),
                    SymbolicUtils.metadata(b),
                )
                @test rebuilt === b
                subbed = Symbolics.substitute(b, Dict(SymbolicUtils.unwrap(_u) => 2.0))
                @test subbed === atom.raw(SymbolicUtils.unwrap(2.0 * _v))
            end

            @testset "lowering round-trips onto the same atom" begin
                @test isequal(c, _to_cnum(to_num(c)))
                @test isequal(_conj_cnum(c), _to_cnum(to_num(_conj_cnum(c))))
                # negative exponents lower to a readable inverse that re-recognizes back
                inv_c = _conj_cnum(c)
                @test isequal(inv_c, _to_cnum(to_num(inv_c)))
            end

            @testset "conjugation law" begin
                @test !isequal(conj(c), c)
                if atom.unimodular
                    @test conj(c) * c == _CNUM_ONE
                    @test _mul_cnum(_mul_cnum(c, c), _conj_cnum(_mul_cnum(c, c))) ==
                        _CNUM_ONE
                end
                # the operator-level adjoint must agree with the coefficient-level one
                @test isequal(_to_cnum(qadjoint(Num(b))), _conj_cnum(c))
            end

            @testset "complex scalars keep it whole" begin
                # `im * c` must stay one factor; splitting it is what loses cancellation.
                @test _is_poly(im * c)
                @test (im * c) * conj(c) == _to_cnum(im)
            end

            @testset "a literal argument folds to a number" begin
                @test _is_native(atom.build(Num(0.75)))
                @test _fold_const(SymbolicUtils.unwrap(to_num(atom.build(Num(0.75))))) isa
                    Complex
            end

            @testset "derivative rule is registered" begin
                # Without one, `expand_derivatives` hits the global `nothing` fallback and
                # leaves an inert `Differential` node that nothing downstream can evaluate.
                d = Symbolics.expand_derivatives(Symbolics.Differential(_u)(Num(b)))
                @test !occursin("Differential", string(d))
            end

            @testset "both renderers know the head" begin
                # `show` and Latexify are separate code paths; a head registered in only one
                # renders as its raw internal name in the other.
                h = FockSpace(:f)
                op = Destroy(h, :a)
                shown = sprint(show, c * op)
                tex = String(latexify(c * op))
                @test !occursin(atom.name, shown)
                @test !occursin(atom.name, tex)
                # a conjugate pair folds the same way in both
                pair = (c + conj(c)) * 0.5
                @test occursin("cos", sprint(show, pair * op)) ==
                    occursin("cos", String(latexify(pair * op)))
            end

            if atom.bounded
                @testset "bounded atoms admit a finite sup-norm" begin
                    # `is_canonical` bounds a floating-point residual instead of sampling it;
                    # an atom with no bound makes every such residual read as nonzero.
                    @test _sup_norm(real(to_num(c))) !== nothing
                end
            end
        end
    end
end
