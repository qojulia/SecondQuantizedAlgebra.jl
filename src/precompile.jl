using PrecompileTools: @setup_workload, @compile_workload

# Targeted precompile directives for internal methods.
precompile(_ordered_qadd, (CNum, Vector{QSym}))
precompile(_site_sort!, (Vector{QSym},))
precompile(_addto!, (QTermDict, Vector{QSym}, CNum))
precompile(_iszero_cnum, (CNum,))
precompile(_neg_cnum, (CNum,))
precompile(_merge_unique, (Vector{Index}, Vector{Index}))
precompile(_merge_unique, (Vector{Tuple{Index, Index}}, Vector{Tuple{Index, Index}}))

# Minimal workload: a single commutator covers Sort, Dict, Symbolics Num paths,
# ordering, QAdd subtraction, and _merge_unique.
@setup_workload begin
    @compile_workload begin
        h = FockSpace(:a)
        a = Destroy(h, :a)
        ad = Create(h, :a)
        commutator(a, ad)
    end
end
# precompile load from 6s to 11s, but TTFX for mul, commutator below 1s.
# do not add simplify as this adds an additional 7s of precompilation time due to SU.jl TTFX.
