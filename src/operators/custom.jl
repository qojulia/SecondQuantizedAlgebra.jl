"""
    OperatorKind(name::Symbol, space::Type{<:HilbertSpace}, n::Int; adjoint, relations, realization, suffixes, latex)

Operator kind `name` with `n` generators acting on spaces of type `space`, for an algebra
that no built-in operator covers.

Generators are numbered `1:n` and made with [`CustomOperator`](@ref). The keywords define
the kind:

  - `adjoint`: the generator `adjoint[g]` is the adjoint of generator `g`.
  - `relations`: a function that receives the `n` generators and returns pairs
    `(x, y) => xy` that rewrite the product `x * y`. Products without a relation stay as
    written. Each right-hand side is a number or an expression in the generators, and
    rewriting must bring every product closer to an ordered form, or it does not terminate.
  - `realization`: `n` expressions in built-in operators on a `space`, with numeric
    coefficients, used by [`to_numeric`](@ref). Without it, conversion throws.
  - `suffixes`, `latex`: what each generator prints after its name, by default its number.

The relations apply to generators with the same name and site; all other operators act
independently of them. Creating a kind with a name in use redefines it. Like operator names,
a kind lives in the current session, so a package creates its kinds in `__init__`.

See also [`CustomOperator`](@ref), [`operator_kind`](@ref).
"""
struct OperatorKind
    id::UInt16
end

struct KindData
    name::Symbol
    space::DataType
    suffixes::Vector{String}
    latex::Vector{String}
    adjoint::Vector{Int32}
    ordered::Matrix{Bool}
    rules::Matrix{Vector{Tuple{CNum, Vector{Int32}}}}
    realization::Vector{Vector{Tuple{ComplexF64, Vector{Op}}}}
    probe::Op
end

const KIND_DATA = KindData[]
const KIND_RANK = Int32[]
const KIND_ID = Dict{Symbol, UInt16}()

kind_data(o::Op) = KIND_DATA[o.custom]
kind_rank(id::UInt16) = id == 0 ? Int32(0) : KIND_RANK[id]

Base.show(io::IO, k::OperatorKind) = print(io, KIND_DATA[k.id].name)

generator(id::UInt16, name::Symbol, g::Int) =
    Op(OP_CUSTOM, id, name_id(name), Int32(1), NO_INDEX, Int32(g), Int32(0), Int32(0), Int32(0))

no_relations(::Vector{Op}) = ()

subscript_string(g::Int) = sprint(write_subscript, g)

function OperatorKind(
        name::Symbol, space::Type{<:HilbertSpace}, n::Int;
        adjoint, relations = no_relations, realization = (),
        suffixes = map(subscript_string, 1:n), latex = map(g -> string("_{", g, "}"), 1:n),
    )
    1 <= n <= typemax(Int32) || throw(ArgumentError("an operator kind needs at least one generator, got $n"))
    adj = Int32[g for g in adjoint]
    (length(adj) == n && all(g -> 1 <= g <= n, adj) && all(g -> adj[adj[g]] == g, 1:n)) ||
        throw(ArgumentError("`adjoint` must map the $n generators onto themselves and be its own inverse, got $adjoint"))
    sfx = String[s for s in suffixes]
    ltx = String[s for s in latex]
    (length(sfx) == n && length(ltx) == n) || throw(ArgumentError("`suffixes` and `latex` need one entry per generator"))
    return lock(INTERN_LOCK) do
        id = get(KIND_ID, name, 0x0000)
        free(id) = KindData(
            name, space, sfx, ltx, adj, fill(true, n, n),
            [Tuple{CNum, Vector{Int32}}[] for _ in 1:n, _ in 1:n],
            Vector{Tuple{ComplexF64, Vector{Op}}}[], generator(id, name, 1),
        )
        if id == 0
            id = UInt16(length(KIND_DATA) + 1)
            push_rank!(name)
            push!(KIND_DATA, free(id))
            KIND_ID[name] = id
        else
            KIND_DATA[id] = free(id)
        end
        gens = Op[generator(id, name, g) for g in 1:n]
        ordered, rules = kind_rules(relations(gens), gens)
        real, probe = kind_realization(realization, n, gens[1])
        KIND_DATA[id] = KindData(name, space, sfx, ltx, adj, ordered, rules, real, probe)
        return OperatorKind(id)
    end
end

function push_rank!(name::Symbol)
    s = string(name)
    r = Int32(1)
    for j in eachindex(KIND_RANK)
        string(KIND_DATA[j].name) < s ? (r += Int32(1)) : (KIND_RANK[j] += Int32(1))
    end
    return push!(KIND_RANK, r)
end

function kind_rules(relations, gens::Vector{Op})
    n = length(gens)
    ordered = fill(true, n, n)
    rules = [Tuple{CNum, Vector{Int32}}[] for _ in 1:n, _ in 1:n]
    for relation in relations
        x, y = relation.first
        (x isa Op && y isa Op && x in gens && y in gens) ||
            throw(ArgumentError("a relation must start from a pair of generators of the kind, got $(relation.first)"))
        ordered[x.l1, y.l1] || throw(ArgumentError("two relations for the product $x * $y"))
        ordered[x.l1, y.l1] = false
        for (term, c) in relation_qadd(relation.second, x).arguments
            push!(rules[x.l1, y.l1], (c, Int32[o.l1 for o in term.ops]))
        end
    end
    return ordered, rules
end

relation_qadd(x::Number, ::Op) = single_qadd(to_cnum(x), Op[])
relation_qadd(x::Op, g::Op) = relation_qadd(to_qadd(x), g)
function relation_qadd(q::QAdd, g::Op)
    isempty(q.indices) || relation_error("a sum over an index")
    for (term, _) in q.arguments
        isempty(term.ne) || relation_error("an index constraint")
        for o in term.ops
            (o.kind === OP_CUSTOM && o.custom == g.custom && same_site(o, g)) ||
                relation_error("the operator $o, which is not a generator of the kind")
        end
    end
    return q
end
relation_qadd(x, ::Op) = relation_error("a $(typeof(x))")

@noinline relation_error(what::String) =
    throw(ArgumentError("the right-hand side of a relation must be a number or an expression in the generators, got $what"))

function kind_realization(realization, n::Int, template::Op)
    real = Vector{Tuple{ComplexF64, Vector{Op}}}[]
    probe = template
    isempty(realization) && return real, probe
    length(realization) == n || throw(ArgumentError("`realization` needs one expression per generator, got $(length(realization))"))
    for x in realization
        q = x isa Number ? single_qadd(to_cnum(x), Op[]) : to_qadd(x)
        isempty(q.indices) || realization_error("a sum over an index")
        terms = Tuple{ComplexF64, Vector{Op}}[]
        for (term, c) in q.arguments
            coeff_is_const(to_num(c)) || realization_error("the symbolic coefficient $(to_num(c))")
            for o in term.ops
                (o.kind !== OP_CUSTOM && !has_index(o.index) && o.space_index == 1) ||
                    realization_error("the operator $o, which is not a built-in operator on a single space")
                probe === template && (probe = o)
            end
            push!(terms, (to_complex(c), copy(term.ops)))
        end
        push!(real, terms)
    end
    probe === template && realization_error("no operator")
    return real, probe
end

@noinline realization_error(what::String) =
    throw(ArgumentError("a realization must be numeric expressions in built-in operators on one space, got $what"))

"""
    CustomOperator(kind::OperatorKind, h::HilbertSpace, name::Symbol, g::Int)
    CustomOperator(kind::OperatorKind, h::ProductSpace, name::Symbol, g::Int, [subsystem::Int])

Generator `g` of the operator kind `kind`, named `name` and acting on `h`.

Generators of one kind with the same name and site follow the relations of the kind;
different names are independent, as two modes are. On a [`ProductSpace`](@ref),
`subsystem` selects the factor and may be omitted when only one factor has the space type
of `kind`. The result is an ordinary [`Op`](@ref) with [`optype`](@ref) `OP_CUSTOM`, so it
can be indexed, summed, averaged and converted like any other operator.

See also [`OperatorKind`](@ref), [`operator_label`](@ref).
"""
function CustomOperator(kind::OperatorKind, h::HilbertSpace, name::Symbol, g::Int)
    data = KIND_DATA[kind.id]
    h isa data.space || throw(ArgumentError("operator kind $kind acts on a $(nameof(data.space)), got a $(nameof(typeof(h)))"))
    1 <= g <= length(data.adjoint) || throw(ArgumentError("operator kind $kind has $(length(data.adjoint)) generators, got $g"))
    return generator(kind.id, name, g)
end
function CustomOperator(kind::OperatorKind, h::ProductSpace, name::Symbol, g::Int, subsystem::Int = kind_factor(kind, h))
    1 <= subsystem <= length(h.spaces) ||
        throw(ArgumentError("Index $subsystem out of range for ProductSpace with $(length(h.spaces)) spaces"))
    return set_acts_on(CustomOperator(kind, h.spaces[subsystem], name, g), subsystem)
end
CustomOperator(::OperatorKind, ::HilbertSpace, name::AbstractString, args...) = name_must_be_symbol(name)

function kind_factor(kind::OperatorKind, h::ProductSpace)
    space = KIND_DATA[kind.id].space
    n = count(s -> s isa space, h.spaces)
    n == 1 || throw(ArgumentError("$n $(nameof(space)) subspaces found in ProductSpace; specify which one with the subspace index argument"))
    return position_of(s -> s isa space, h.spaces)
end

is_custom(o::Op) = o.kind === OP_CUSTOM

custom_adjoint(o::Op) = with_levels(o, kind_data(o).adjoint[o.l1], 0)

"""
    operator_kind(o::Op) -> OperatorKind

The kind of an operator made with [`CustomOperator`](@ref). Throws an `ArgumentError` for
any other operator.

See also [`operator_label`](@ref), [`optype`](@ref).
"""
function operator_kind(o::Op)
    is_custom(o) || wrong_kind(:operator_kind, "a CustomOperator", o)
    return OperatorKind(o.custom)
end

"""
    operator_label(o::Op) -> Int

The generator number of an operator made with [`CustomOperator`](@ref). Throws an
`ArgumentError` for any other operator.

See also [`operator_kind`](@ref).
"""
function operator_label(o::Op)
    is_custom(o) || wrong_kind(:operator_label, "a CustomOperator", o)
    return Int(o.l1)
end
