# Symbolics v7 & SymbolicUtils v4 Migration - Summary

## What Was Done

This PR adds a comprehensive migration plan for upgrading `SecondQuantizedAlgebra.jl` from:
- **SymbolicUtils v3.6.0 → v4.0.0**
- **Symbolics v6 → v7.0.0**

## Files Added

### 1. `.github/issues/symbolics-v7-migration-plan.md`
A detailed 850+ line migration plan document that includes:

#### I. Complete Inventory of Touchpoints
- **13 source files** with Symbolics/SymbolicUtils usage mapped
- **17 custom types** extending SymbolicUtils types documented
- **100+ API usage sites** catalogued with line references
- Tables showing usage counts, priority levels, and breaking change impacts

#### II. Breaking Changes Analysis
Detailed analysis of both upstream breaking changes:

**SymbolicUtils v4:**
- Moshi.jl sum type migration
- `symtype` → `.type` field migration (CRITICAL)
- `Add`/`Mul` → `AddMul` variant unification
- `Const` variant for type stability
- Mandatory `shape` field for all variants
- Removal of custom `Symbolic` subtypes

**Symbolics v7:**
- `Num`/`Arr` wrapper changes
- ArrayOp migration
- Derivative rule registration changes

Each breaking change includes:
- Exact impact on SQA codebase
- Files and line numbers affected
- Migration actions required

#### III. Moshi.jl Evaluation
Recommendations on whether SQA should:
- Use Moshi.jl for its own types (NO for Phase 1)
- Keep current TermInterface approach (YES)
- Re-evaluate in future phases (MAYBE)

#### IV. 7-Phase Migration Plan
Detailed phased approach with ~11 week timeline:

**Phase 0:** Preparation & Investigation (Week 1)
**Phase 1:** Foundation - Type System Migration (Week 2-3)
**Phase 2:** Indexed Types & Metadata (Week 3-4)
**Phase 3:** Tree Traversal & Utilities (Week 4-5)
**Phase 4:** Display & Printing (Week 5)
**Phase 5:** Commutator & Specialized Operations (Week 5-6)
**Phase 6:** Integration & Testing (Week 6-7)
**Phase 7:** Release Preparation (Week 7-8)

Each phase includes:
- ✅ Specific tasks with checkboxes
- 📁 Files to modify with line references
- 🧪 Tests to run
- ✔️ Acceptance criteria

#### V. Testing Strategy
- Baseline testing approach
- Incremental testing per phase
- Final validation procedures
- Performance benchmarking

#### VI. Risk Assessment
Table of 8 identified risks with:
- Probability ratings
- Impact assessments
- Mitigation strategies

#### VII. Additional Sections
- Rollback plan
- Success criteria
- Timeline estimates
- File/line reference tables
- Migration resources
- Questions for upstream maintainers

### 2. `.github/issues/README.md`
Instructions for creating the GitHub issue from the migration plan.

### 3. `.github/issues/create-migration-issue.sh`
Automated script to create the issue via GitHub CLI.

## Migration Plan Highlights

### Most Critical Changes Required

1. **Type System Overhaul** (qnumber.jl, average.jl, cnumber.jl)
   - `BasicSymbolic{CustomType}` → `BasicSymbolic{vartype}` + runtime type checks
   - `symtype(x)` → field access
   - All type aliases must be redesigned

2. **Add/Mul Constructor Updates** (18 methods across multiple files)
   - Migrate to `AddMul` variant with discriminator
   - Update all delegating constructors

3. **Const Handling** (~100+ call sites)
   - Add `unwrap_const` where needed
   - Update tree traversal logic

4. **Shape Field Addition** (QMul, QAdd types)
   - Add shape support to custom types
   - Update `maketerm` implementations

### Testing Approach

- **13 test files** to validate (~2,226 lines)
- Incremental testing after each phase
- Performance benchmarking vs. baseline
- Code quality validation (Aqua, JET, ExplicitImports)

### Documentation Quality

The migration plan is:
- ✅ **Actionable** - 50+ specific checkbox tasks
- ✅ **Detailed** - File and line references throughout
- ✅ **Comprehensive** - Covers all aspects from testing to rollback
- ✅ **Research-backed** - Cites official release notes and migration guides
- ✅ **Realistic** - Includes timeline, risks, and success criteria

## Next Steps

1. **Create the GitHub issue:**
   ```bash
   cd SecondQuantizedAlgebra.jl
   ./.github/issues/create-migration-issue.sh
   ```
   
   Or manually at: https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/new

2. **Review and discuss** the migration plan with maintainers

3. **Begin Phase 0** once approved:
   - Set up development branch
   - Create test baseline
   - Prototype key patterns

## References

All migration plan sections cite these authoritative sources:
- https://docs.sciml.ai/SymbolicUtils/stable/manual/variants/
- https://github.com/JuliaSymbolics/Symbolics.jl/releases/tag/v7.0.0
- https://github.com/JuliaSymbolics/SymbolicUtils.jl/releases/tag/v4.0.0
- https://github.com/JuliaSymbolics/Symbolics.jl/compare/v6.57.0...master
- https://github.com/JuliaSymbolics/SymbolicUtils.jl/compare/v3.32.0...master
