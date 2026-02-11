# Migration Plan: Upgrade to Symbolics v7 & SymbolicUtils v4

## Overview

This issue tracks the migration of `SecondQuantizedAlgebra.jl` to:
- **SymbolicUtils v4** (currently v3.6.0)
- **Symbolics v7** (currently v6)

Both upgrades involve **breaking changes** to the symbolic algebra system, particularly around type representations, variant structures, and internal APIs. This migration will require careful refactoring across the entire codebase.

## Migration Resources

- [SymbolicUtils v4 Variants Documentation](https://docs.sciml.ai/SymbolicUtils/stable/manual/variants/)
- [Symbolics v7.0.0 Release Notes](https://github.com/JuliaSymbolics/Symbolics.jl/releases/tag/v7.0.0)
- [SymbolicUtils v4.0.0 Release Notes](https://github.com/JuliaSymbolics/SymbolicUtils.jl/releases/tag/v4.0.0)
- [Symbolics v6.57.0...master diff](https://github.com/JuliaSymbolics/Symbolics.jl/compare/v6.57.0...master)
- [SymbolicUtils v3.32.0...master diff](https://github.com/JuliaSymbolics/SymbolicUtils.jl/compare/v3.32.0...master)

---

## I. Current Symbolics/SymbolicUtils Touchpoints

### A. Core Files with Heavy Symbolic Usage

| File | Lines | Key APIs Used | Priority |
|------|-------|---------------|----------|
| `src/qnumber.jl` | 260 | `symtype`, `promote_symtype`, `iscall`, `operation`, `arguments`, `TermInterface.*`, `hash` | **CRITICAL** |
| `src/average.jl` | 100 | `BasicSymbolic`, `Term{T}`, `Sym{FnType}`, `promote_symtype`, `symtype`, `Add`, `Mul` | **CRITICAL** |
| `src/cnumber.jl` | 180 | `Sym{Complex{Real}}`, `Sym{Real}`, `setmetadata`, `Symbolic{T}`, `BasicSymbolic{T}` | **CRITICAL** |
| `src/index_average.jl` | 350 | `BasicSymbolic{T}`, `Sym{T}`, `setmetadata`, `hasmetadata`, `metadata`, `Add`, `Mul`, `_iszero` | **HIGH** |
| `src/indexing.jl` | 380 | `Sym{T}`, `setmetadata`, `BasicSymbolic`, `Mul`, `Div`, `iscall`, `operation`, `arguments` | **HIGH** |
| `src/index_double_sums.jl` | 130 | `iscall`, `arguments`, `operation` | **MEDIUM** |
| `src/index_utils.jl` | 100 | `BasicSymbolic`, `iscall`, `operation`, `arguments`, `maketerm`, `metadata` | **MEDIUM** |
| `src/index_numbered_operator.jl` | 35 | `iscall`, `operation`, `arguments` | **MEDIUM** |
| `src/utils.jl` | 200 | `operation`, `arguments`, `maketerm`, `metadata`, `_iszero`, `_isone` | **MEDIUM** |
| `src/printing.jl` | 60 | `show_term`, `arguments`, `operation` | **LOW** |
| `src/latexify_recipes.jl` | 140 | `iscall`, `operation`, `arguments` | **LOW** |
| `src/commutator.jl` | 200 | Direct `.arguments` field access | **LOW** |

**Total touchpoints:** 13 files, ~2200+ lines of symbolic code

### B. Custom Types Extending SymbolicUtils

#### 1. Types Using `Sym{T}` Construction Pattern

| Type | File | Lines | Constructor Pattern | Breaking Change Impact |
|------|------|-------|---------------------|------------------------|
| `Parameter` | cnumber.jl | L14-19 | `Sym{Complex{Real}}` + `setmetadata` | **HIGH** - symtype → vartype + type field |
| `RealParameter` | cnumber.jl | L118-123 | `Sym{Real}` + `setmetadata` | **HIGH** - symtype → vartype + type field |
| `AvgSym` & `sym_average` | average.jl | L13-16 | `Sym{FnType{Tuple{QNumber},AvgSym}}` | **CRITICAL** - FnType compatibility unknown |
| `IndexedVariable` | indexing.jl | L49-58 | `Sym{IndexedVariable}` + `setmetadata` | **HIGH** - recursive symtype |
| `DoubleIndexedVariable` | indexing.jl | L75-86 | `Sym{DoubleIndexedVariable}` + `setmetadata` | **HIGH** - recursive symtype |
| `IndexedAverageSum` | index_average.jl | L40-45 | `Sym{IndexedAverageSum}` + double `setmetadata` | **HIGH** - complex metadata structure |
| `IndexedAverageDoubleSum` | index_average.jl | L89-99 | `Sym{IndexedAverageDoubleSum}` + `setmetadata` | **HIGH** |
| `SpecialIndexedAverage` | index_average.jl | L213-214 | `Sym{SpecialIndexedAverage}` + `setmetadata` | **HIGH** |

#### 2. Types Using `BasicSymbolic{T}` Type Aliases

| Type Alias | File | Definition | Usage Count | Breaking Change Impact |
|------------|------|------------|-------------|------------------------|
| `Average` | average.jl | L11: `BasicSymbolic{<:AvgSym}` | 40+ | **CRITICAL** - vartype affects generic parameter |
| `symbolics_terms` | index_average.jl | `Union{<:Average,<:BasicSymbolic{<:CNumber}}` | 10+ | **HIGH** |
| Various dispatches | cnumber.jl | L157-174: Methods on `BasicSymbolic{Complex{RNumber}}`, `BasicSymbolic{CNumber}` | 6+ | **MEDIUM** |

#### 3. Custom QNumber Type Hierarchy (TermInterface-based)

| Type | Parent | TermInterface Impl | Breaking Change Impact |
|------|--------|-------------------|------------------------|
| `QNumber` | abstract | `head`, `metadata` | **LOW** - TermInterface stable |
| `QSym` | `QNumber` | `iscall=false`, custom `hash` | **LOW** |
| `QTerm` | `QNumber` | `iscall=true`, `operation`, `arguments`, `maketerm` | **MEDIUM** - shape field required |
| `QMul{M}` | `QTerm` | L114-122: Custom `maketerm`, `metadata` | **MEDIUM** |
| `QAdd` | `QTerm` | L238-241: Custom `maketerm`, `NO_METADATA` | **MEDIUM** |

### C. SymbolicUtils API Usages by Category

#### Type System APIs (All BREAKING in v4)

| API | Current Usage | v4 Change | Files Affected |
|-----|---------------|-----------|----------------|
| `symtype(x)` | Returns type of symbolic value | **REMOVED** - Now stored in `.type` field, not type parameter | qnumber.jl L59, average.jl L22 |
| `promote_symtype(f, ...)` | Type promotion for operations | **MODIFIED** - Must account for new vartype system | qnumber.jl L51-57, average.jl L19 |
| `BasicSymbolic{T}` | Generic parameter `T` is symtype | **CHANGED** - `T` now represents vartype (SymReal/SafeReal/TreeReal) | average.jl L11, index_average.jl (many), cnumber.jl L157+ |

#### Variant Construction APIs (BREAKING in v4)

| API | Current Usage | v4 Change | Files Affected |
|-----|---------------|-----------|----------------|
| `Term{T}(op, args)` | Construct function application terms | **COMPATIBLE** - Still exists but `T` is vartype, not symtype | average.jl L26 |
| `Sym{T}(name)` | Create symbolic variables | **MODIFIED** - `T` interpretation changed; needs vartype consideration | cnumber.jl L16/L119, indexing.jl L53/L81, index_average.jl (many) |
| `Add(Type{T}, coeff, dict; kw...)` | Create addition expressions | **MERGED** into `AddMul` variant with `AddMulVariant` discriminator | average.jl L43, index_average.jl (8 methods) |
| `Mul(Type{T}, coeff, dict; kw...)` | Create multiplication expressions | **MERGED** into `AddMul` variant with `AddMulVariant` discriminator | average.jl L46, index_average.jl (8 methods) |

#### Tree Traversal APIs (MOSTLY STABLE)

| API | Usage Count | v4 Status | Notes |
|-----|-------------|-----------|-------|
| `iscall(x)` | 30+ | ✅ STABLE | No changes |
| `operation(x)` | 40+ | ✅ STABLE | No changes |
| `arguments(x)` | 100+ | ⚠️ MODIFIED | Now always returns `BasicSymbolic` (wrapped in `Const` if not symbolic) |
| `TermInterface.maketerm` | 10+ | ✅ STABLE | No changes (TermInterface is separate) |
| `TermInterface.metadata` | 20+ | ⚠️ POSSIBLY MODIFIED | Shape field now mandatory in all variants |

#### Utility APIs (MODIFIED in v4)

| API | Usage Count | v4 Change | Impact |
|-----|-------------|-----------|--------|
| `_iszero(x)` | 15+ | ⚠️ May need `unwrap_const` for Const handling | index_average.jl, qnumber.jl, utils.jl |
| `_isone(x)` | 10+ | ⚠️ May need `unwrap_const` for Const handling | qnumber.jl, utils.jl |
| `setmetadata(x, k, v)` | 12+ | ❓ UNKNOWN - Check if shape compatibility needed | All Sym{T} constructors |
| `hasmetadata(x, k)` | 3 | ✅ LIKELY STABLE | index_average.jl L268-272 |

#### Hash/Equality (MODIFIED in v4)

| API | Usage | v4 Change |
|-----|-------|-----------|
| `hash(::QAdd, ...)` | qnumber.jl L229 | Must call `SymbolicUtils.hashvec` - verify v4 availability |
| `hash(::QMul, ...)` | qnumber.jl L109 | Must call `SymbolicUtils.hashvec` - verify v4 availability |

---

## II. Expected Breaking Changes

### A. SymbolicUtils v4 Breaking Changes (from release notes)

#### 1. **Moshi.jl Sum Type Migration** ⚠️ CRITICAL

**Change:** `BasicSymbolic` now uses [Moshi.jl](https://github.com/JuliaComputing/Moshi.jl) for variant representation instead of Unityper.

**Impact on SQA:**
- Pattern matching on variants requires `Moshi.Match.@match` instead of previous methods
- All `BasicSymbolic` handling code must be reviewed
- Possible performance implications (positive or negative)

**Files affected:** All files using `BasicSymbolic{T}` dispatch

**Migration action required:**
- [ ] Evaluate if custom pattern matching exists that needs `@match` conversion
- [ ] Test performance of Moshi-based variants vs. current Unityper

#### 2. **symtype → type Field Migration** 🔴 BREAKING - CRITICAL

**Change:** 
- `BasicSymbolic{T}` parameter `T` no longer represents symtype
- `T` now represents **vartype** (algebra type: `SymReal`, `SafeReal`, `TreeReal`)
- Actual symtype stored in `.type` field of each variant

**Impact on SQA:**
- ❌ **BROKEN:** `const Average = BasicSymbolic{<:AvgSym}` - `AvgSym` is not a valid vartype
- ❌ **BROKEN:** All `BasicSymbolic{Complex{RNumber}}` dispatches
- ❌ **BROKEN:** All `BasicSymbolic{CNumber}` dispatches
- ❌ **BROKEN:** `symtype(x::T) where {T<:QNumber} = T` method (L59 in qnumber.jl)
- ❌ **BROKEN:** `promote_symtype` methods may need redesign

**Files affected (CRITICAL):**
- `src/average.jl` L11, L22 - Average type alias and symtype method
- `src/qnumber.jl` L59 - symtype implementation
- `src/cnumber.jl` L157-174 - BasicSymbolic method dispatches
- `src/index_average.jl` - All `BasicSymbolic{IndexedAverageSum}` usages (~20 instances)
- `src/index_utils.jl` L84-91 - Generic `BasicSymbolic` reconstruction

**Migration action required:**
- [ ] Redefine type aliases to use vartype + separate type checking
- [ ] Replace `symtype(x)` calls with field access: `SymbolicUtils.gettype(x)` (or equivalent v4 API)
- [ ] Update all `BasicSymbolic{CustomType}` dispatches to use vartype + runtime type checks
- [ ] Verify `promote_symtype` still makes sense or needs replacement

#### 3. **Add/Mul → AddMul Variant Unification** 🔴 BREAKING - HIGH

**Change:** `Add` and `Mul` merged into single `AddMul` variant with `AddMulVariant.T` enum field (from EnumX.jl).

**Impact on SQA:**
- ❌ **BROKEN:** 8 methods calling `SymbolicUtils.Add(::Type{CustomType}, ...)` in `index_average.jl`
- ❌ **BROKEN:** 8 methods calling `SymbolicUtils.Mul(::Type{CustomType}, ...)` in `index_average.jl`
- ❌ **BROKEN:** 2 methods in `average.jl` (L43-47)
- ⚠️ **UNKNOWN:** Will `Add(CNumber, ...)` still work as constructor, or must we use `AddMul(..., variant=AddMulVariant.Add)`?

**Files affected:**
- `src/average.jl` L43-47
- `src/index_average.jl` L328-367 (16 methods redirecting Add/Mul to CNumber symtype)

**Migration action required:**
- [ ] Determine v4 API for constructing Add vs. Mul (likely `AddMul` with variant parameter)
- [ ] Update all 18 Add/Mul constructor calls
- [ ] Verify that delegation pattern (redirect to CNumber) still makes sense

#### 4. **Const Variant for Type Stability** ⚠️ MEDIUM

**Change:** Non-symbolic constants wrapped in `Const` variant for type stability.

**Impact on SQA:**
- `arguments(x)` now returns `BasicSymbolic[]` even for constant arguments (wrapped in `Const`)
- Must use `unwrap_const(arg)` to extract actual value
- Affects all tree traversal and argument processing

**Files affected:**
- All files using `arguments(...)`: qnumber.jl, average.jl, index_utils.jl, utils.jl, printing.jl, latexify_recipes.jl, commutator.jl, index_average.jl, indexing.jl, index_double_sums.jl, index_numbered_operator.jl

**Migration action required:**
- [ ] Audit all `arguments(x)` usages (~100+ call sites)
- [ ] Add `unwrap_const` where needed for constant extraction
- [ ] Update `_iszero` / `_isone` checks to handle `Const` wrapping

#### 5. **Mandatory shape Field** 🔴 BREAKING - MEDIUM

**Change:** All variants now have a `shape` field for first-class symbolic array support.

**Impact on SQA:**
- Custom types like `QMul`, `QAdd` must provide shape information
- `TermInterface.maketerm` calls may need shape parameter
- Metadata handling may interact with shape

**Files affected:**
- `src/qnumber.jl` L114-120, L240 - QMul/QAdd maketerm implementations

**Migration action required:**
- [ ] Determine default shape for QNumber types (likely scalar: `SymbolicUtils.NotSymbolic()`)
- [ ] Update `maketerm` implementations to specify shape
- [ ] Test that shape propagation doesn't break existing behavior

#### 6. **Removal of `Symbolic` Custom Subtype Support** ⚠️ LOW

**Change:** SymbolicUtils no longer allows defining custom `Symbolic` subtypes.

**Impact on SQA:**
- ✅ **OK** - SQA doesn't define custom `Symbolic` subtypes
- Only uses `Symbolic{T}` for dispatch in a few places (cnumber.jl L28, L132)

**Files affected:**
- `src/cnumber.jl` L28, L132 - `Base.adjoint(x::SymbolicUtils.Symbolic{<:CNumber})`

**Migration action required:**
- [ ] Verify that dispatch on `Symbolic{T}` still works for standard variants
- [ ] If removed, replace with `BasicSymbolic` dispatch + type field check

#### 7. **Removal of PolyForm** ✅ NONE

**Change:** `PolyForm` removed; use `to_poly!` / `from_poly` instead.

**Impact on SQA:** ✅ No usage of PolyForm detected in codebase.

### B. Symbolics v7 Breaking Changes (from release notes)

#### 1. **Num/Arr Wrapper Changes** ⚠️ LOW-MEDIUM

**Change:**
- `Num` and `Arr` now wrap `BasicSymbolic{SymReal}` (instead of `::Any`)
- `unwrap(x)` returns `BasicSymbolic` (not raw value)
- New function `Symbolics.value(x)` = `unwrap_const(unwrap(x))` recovers old `unwrap` behavior

**Impact on SQA:**
- ✅ **LIKELY OK** - SQA doesn't appear to use `Num`/`Arr` directly
- No `unwrap` calls detected in source
- Primarily works with raw `BasicSymbolic` and custom types

**Migration action required:**
- [ ] Verify no hidden `Num`/`Arr` usage in tests or dependencies
- [ ] If `unwrap` is added in future, use `Symbolics.value` instead

#### 2. **ArrayOp Migration to SymbolicUtils** ✅ NONE

**Change:** ArrayOp moved to SymbolicUtils; ArrayMaker temporarily removed.

**Impact on SQA:** ✅ No usage of ArrayOp or ArrayMaker detected.

#### 3. **@register_array_symbolic ndims Requirement** ⚠️ LOW

**Change:** `@register_array_symbolic` now soft-requires `ndims` specification.

**Impact on SQA:** ✅ No usage of `@register_array_symbolic` detected.

#### 4. **Derivative Rule Registration Changes** ✅ NONE

**Change:** `Symbolics.derivative` methods → `@register_derivative` macro.

**Impact on SQA:** ✅ No custom derivative rules detected.

#### 5. **CallWithMetadata Removal** ⚠️ LOW

**Change:** `Symbolics.CallWithMetadata` removed; callable symbolics auto-transfer metadata. `CallAndWrap` added.

**Impact on SQA:** ✅ No usage of `CallWithMetadata` detected. May benefit from auto-metadata-transfer.

---

## III. Moshi.jl Consideration

### What is Moshi.jl?

Moshi.jl is a **sum type** library for Julia providing:
- Efficient tagged union types (similar to Rust enums)
- Pattern matching via `@match` macro
- Type-stable variant representation

SymbolicUtils v4 uses Moshi internally for the `BasicSymbolic` sum type (replacing Unityper).

### Current Status in SQA

**Direct usage:** ✅ None detected
**Indirect usage:** ⚠️ Via SymbolicUtils v4 (if migrated)

### Should SQA Use Moshi Directly?

#### Option A: Use Moshi for Custom Types (e.g., QTerm Variants)

**Pros:**
- Type stability for QMul, QAdd, SingleSum, DoubleSum variants
- Pattern matching via `@match` for tree traversal
- Consistency with SymbolicUtils v4 architecture
- Potential performance improvements

**Cons:**
- Adds new dependency
- Learning curve for contributors
- Current TermInterface-based approach works well
- May not be necessary for relatively simple type hierarchy

#### Option B: Keep Current TermInterface Approach

**Pros:**
- ✅ Minimal changes to existing codebase
- ✅ TermInterface is stable and well-understood
- ✅ QNumber hierarchy is already well-designed
- ✅ No new dependencies

**Cons:**
- May miss performance optimizations
- Less alignment with SymbolicUtils v4 patterns

### **Recommendation**

**Phase 1 (v0.5):** ❌ Do NOT add Moshi for SQA's own types
- Current TermInterface approach is sufficient
- Focus migration effort on SymbolicUtils v4 compatibility
- Avoid unnecessary complexity during breaking migration

**Phase 2+ (Future):** ❓ Re-evaluate after migration complete
- Benchmark QNumber type dispatch performance
- If bottlenecks found, consider Moshi for QTerm variants
- Wait for broader ecosystem adoption

---

## IV. Proposed Migration Plan

### Phase 0: Preparation & Investigation (Week 1)

#### Tasks

- [ ] **0.1** Set up development branch: `feature/symbolics-v7-migration`
- [ ] **0.2** Create comprehensive test baseline on current v3/v6
  - Run full test suite, record results
  - Identify any pre-existing failures (exclude from migration scope)
  - Document expected behavior for critical features
- [ ] **0.3** Set up v4/v7 experimental environment
  - Create isolated Julia environment with SymbolicUtils@4, Symbolics@7
  - Install minimal SQA dependencies
  - Document any immediate compatibility issues
- [ ] **0.4** Deep-dive into v4 APIs
  - Read https://docs.sciml.ai/SymbolicUtils/stable/manual/variants/ thoroughly
  - Identify exact replacement APIs for all breaking changes
  - Create API translation table (v3 → v4)
- [ ] **0.5** Prototype key patterns
  - Test `BasicSymbolic{T}` with vartype in REPL
  - Test `Sym{CustomType}` construction with new semantics
  - Test `AddMul` variant construction
  - Test `Const` handling in tree traversal

**Acceptance Criteria:** 
- Clear understanding of all required API changes
- Prototypes validate feasibility
- No blockers identified

---

### Phase 1: Foundation - Type System Migration (Week 2-3)

**Goal:** Update core type representations and symtype/vartype handling.

#### 1.1 Update `cnumber.jl` - Parameters & CNumber Types

**File:** `src/cnumber.jl`

**Changes:**

- [ ] **1.1.1** Determine appropriate vartype for CNumber types
  - Options: `SymReal` (default), `SafeReal`, or custom?
  - Recommendation: Start with `SymReal`
- [ ] **1.1.2** Update `Parameter` constructor (L15-19)
  - Keep `Sym{Complex{Real}}` or change to `Sym{vartype}` with `.type` field?
  - Ensure `setmetadata` compatible with v4 shape requirements
- [ ] **1.1.3** Update `RealParameter` constructor (L119-123)
  - Same as above for `Sym{Real}`
- [ ] **1.1.4** Fix `BasicSymbolic{CNumber}` dispatches (L157-174)
  - Replace with `BasicSymbolic{vartype}` + runtime type checks
  - OR use TermInterface-based dispatch if type agnostic
- [ ] **1.1.5** Update `promote_rule` if needed (L24, L128)
- [ ] **1.1.6** Test parameter creation and arithmetic

**Files to modify:**
- `/src/cnumber.jl`

**Tests to run:**
- `test/test_cnumber.jl`

**Acceptance Criteria:**
- Parameters construct without errors
- Arithmetic on parameters works
- Metadata preserved
- All test_cnumber.jl tests pass

---

#### 1.2 Update `qnumber.jl` - QNumber Type System

**File:** `src/qnumber.jl`

**Changes:**

- [ ] **1.2.1** Update `promote_symtype` methods (L51-57)
  - Determine if promote_symtype still used in v4 for custom types
  - If yes, adapt to vartype system
  - If no, replace with alternative promotion mechanism
- [ ] **1.2.2** Update `symtype(x::T) where {T<:QNumber}` (L59)
  - Replace `symtype` function if it's removed in v4
  - Implement accessor for `.type` field if needed
- [ ] **1.2.3** Add shape support to QMul (L94+)
  - Add `shape` field to struct or provide via TermInterface
  - Default: `SymbolicUtils.NotSymbolic()` (scalar)
- [ ] **1.2.4** Update `maketerm(::Type{<:QMul}, ...)` (L114-120)
  - Ensure shape parameter handled correctly
- [ ] **1.2.5** Add shape support to QAdd (L225+)
  - Same as QMul
- [ ] **1.2.6** Update `maketerm(::Type{<:QAdd}, ...)` (L240)
- [ ] **1.2.7** Verify hash implementations (L109, L229)
  - Check `SymbolicUtils.hashvec` still exists in v4
- [ ] **1.2.8** Test QNumber operations (multiplication, addition)

**Files to modify:**
- `/src/qnumber.jl`

**Tests to run:**
- `test/test_fock.jl`, `test/test_nlevel.jl`, `test/test_spin.jl`

**Acceptance Criteria:**
- QMul and QAdd construct correctly
- Operations preserve shape
- Hash and equality work
- Operator tests pass

---

#### 1.3 Update `average.jl` - Average Type & Construction

**File:** `src/average.jl`

**Changes:**

- [ ] **1.3.1** Redefine `Average` type alias (L11)
  - OLD: `BasicSymbolic{<:AvgSym}`
  - NEW: `BasicSymbolic{vartype}` with runtime type checking
  - Consider helper: `is_average(x) = gettype(x) == AvgSym`
- [ ] **1.3.2** Update `sym_average` construction (L13-16)
  - Verify `FnType` compatibility with v4
  - May need to construct differently or use Term directly
- [ ] **1.3.3** Update `promote_symtype(::typeof(sym_average), ...)` (L19)
  - Adapt to v4 type system
- [ ] **1.3.4** Update `symtype(::T) where {T<:Average}` (L22)
  - Replace with v4 equivalent (likely field access)
- [ ] **1.3.5** Update `_average` function (L25-27)
  - Verify `Term{AvgSym}` → `Term{vartype}` with `.type = AvgSym`
- [ ] **1.3.6** Fix `Add(::Type{AvgSym}, ...)` (L43-44)
  - Replace with `AddMul` construction
  - Determine correct variant parameter (AddMulVariant.Add)
- [ ] **1.3.7** Fix `Mul(::Type{AvgSym}, ...)` (L46-47)
  - Replace with `AddMul` construction (AddMulVariant.Mul)
- [ ] **1.3.8** Test average construction and simplification

**Files to modify:**
- `/src/average.jl`

**Tests to run:**
- `test/test_average.jl`

**Acceptance Criteria:**
- Averages construct correctly
- Type checking works
- Addition/multiplication of averages works
- Simplification works
- test_average.jl passes

---

### Phase 2: Indexed Types & Metadata (Week 3-4)

**Goal:** Migrate indexed variable types and complex metadata usages.

#### 2.1 Update `indexing.jl` - IndexedVariable & DoubleIndexedVariable

**File:** `src/indexing.jl`

**Changes:**

- [ ] **2.1.1** Update `IndexedVariable` constructor (L49-58)
  - Adapt `Sym{IndexedVariable}` to v4
  - Ensure metadata with shape compatibility
- [ ] **2.1.2** Update `DoubleIndexedVariable` constructor (L75-86)
  - Same as above
- [ ] **2.1.3** Update `Ranges` type union (L4-6)
  - Verify `Mul`, `Div` still valid types
  - May need `AddMul` instead
- [ ] **2.1.4** Update `change_index` function (L240-286)
  - Audit `iscall`, `operation`, `arguments` usage
  - Add `unwrap_const` where needed for constants
- [ ] **2.1.5** Update `arguments(a::SingleSum)` (L345)
  - Handle `Const` wrapping
- [ ] **2.1.6** Test indexed variables

**Files to modify:**
- `/src/indexing.jl`

**Tests to run:**
- `test/test_index_basic.jl`

**Acceptance Criteria:**
- IndexedVariable constructs correctly
- Indexing operations work
- test_index_basic.jl passes

---

#### 2.2 Update `index_average.jl` - IndexedAverageSum & Related Types

**File:** `src/index_average.jl`

**Changes:**

- [ ] **2.2.1** Update `symbolics_terms` type union (L6)
  - Replace `BasicSymbolic{<:CNumber}` with vartype-based check
- [ ] **2.2.2** Update `IndexedAverageSum` construction (L40-45)
  - Ensure double `setmetadata` compatible with v4
  - Add shape if required
- [ ] **2.2.3** Update `IndexedAverageDoubleSum.innerSum` field (L89)
  - Type should be `BasicSymbolic{vartype}`
- [ ] **2.2.4** Update `SpecialIndexedAverage` construction (L213-214)
- [ ] **2.2.5** Fix all `Add(::Type{T}, ...)` methods (L328, 336, 344, 352)
  - 8 total Add methods → AddMul
- [ ] **2.2.6** Fix all `Mul(::Type{T}, ...)` methods (L330, 338, 346, 354)
  - 8 total Mul methods → AddMul
- [ ] **2.2.7** Update `BasicSymbolic{T}` dispatches throughout
  - L87: `undo_average(a::BasicSymbolic{IndexedAverageSum})`
  - L160: `undo_average(a::BasicSymbolic{IndexedAverageDoubleSum})`
  - L190: `undo_average(a::BasicSymbolic{SpecialIndexedAverage})`
  - L243: `get_indices(a::BasicSymbolic{...})`
  - L246: `get_indices(a::BasicSymbolic{...})`
  - L312: `iscall(a::BasicSymbolic{SpecialIndexedAverage})`
  - L315: `iscall(a::BasicSymbolic{IndexedAverageDoubleSum})`
  - Plus equality/arguments methods
  - Replace all with vartype + type checks
- [ ] **2.2.8** Update `hasmetadata` usage (L268-272)
  - Verify API compatibility
- [ ] **2.2.9** Update `_iszero` / `_isone` calls
  - Add `unwrap_const` if needed (L208, L306, L307)
- [ ] **2.2.10** Test indexed average sums

**Files to modify:**
- `/src/index_average.jl`

**Tests to run:**
- `test/test_average_sums.jl`

**Acceptance Criteria:**
- IndexedAverageSum constructs correctly
- Metadata preserved
- Add/Mul operations work
- test_average_sums.jl passes

---

#### 2.3 Update `index_double_sums.jl` & `index_numbered_operator.jl`

**Files:** `src/index_double_sums.jl`, `src/index_numbered_operator.jl`

**Changes:**

- [ ] **2.3.1** Update `DoubleSum` (index_double_sums.jl)
  - Add shape field if needed
  - Update `arguments` method (L128) for Const handling
- [ ] **2.3.2** Update `NumberedOperator` tree traversal (index_numbered_operator.jl L23-31)
  - Audit `iscall`, `operation`, `arguments` usage
- [ ] **2.3.3** Test double sums

**Files to modify:**
- `/src/index_double_sums.jl`
- `/src/index_numbered_operator.jl`

**Tests to run:**
- `test/test_double_sums.jl`

**Acceptance Criteria:**
- DoubleSums construct correctly
- NumberedOperator works
- test_double_sums.jl passes

---

### Phase 3: Tree Traversal & Utilities (Week 4-5)

**Goal:** Update all tree traversal, reconstruction, and utility functions.

#### 3.1 Update `index_utils.jl` - Tree Utilities

**File:** `src/index_utils.jl`

**Changes:**

- [ ] **3.1.1** Update `inorder!(v::BasicSymbolic)` (L84-91)
  - Change to `inorder!(v::BasicSymbolic{vartype})`
  - Ensure type-agnostic or add type checks
- [ ] **3.1.2** Update `inorder!(v::Average)` (L94-100)
  - Adapt Average type check for v4
- [ ] **3.1.3** Update `get_indices` functions (L8-16)
  - Handle Const in arguments
- [ ] **3.1.4** Update `NumberedOperator(x::QTerm, ...)` (L23-28)
  - Handle Const in arguments
- [ ] **3.1.5** Test tree utilities

**Files to modify:**
- `/src/index_utils.jl`

**Tests to run:**
- Related index tests

**Acceptance Criteria:**
- Tree traversal works
- Reconstruction preserves structure
- No errors in index operations

---

#### 3.2 Update `utils.jl` - General Utilities

**File:** `src/utils.jl`

**Changes:**

- [ ] **3.2.1** Update `embed` functions (L89-96)
  - Verify fieldnames approach still works
  - Update for shape if needed
- [ ] **3.2.2** Update tree reconstruction helpers
  - Add unwrap_const where needed
  - Verify maketerm calls
- [ ] **3.2.3** Update conjugation functions (L162-180)
  - Handle Const in arguments
- [ ] **3.2.4** Audit `_iszero` / `_isone` usages
- [ ] **3.2.5** Test utilities

**Files to modify:**
- `/src/utils.jl`

**Tests to run:**
- `test/test_utils.jl`

**Acceptance Criteria:**
- Embedding works
- Conjugation works
- test_utils.jl passes

---

### Phase 4: Display & Printing (Week 5)

**Goal:** Update all display, printing, and latexification.

#### 4.1 Update `printing.jl` & `latexify_recipes.jl`

**Files:** `src/printing.jl`, `src/latexify_recipes.jl`

**Changes:**

- [ ] **4.1.1** Update `show_term(io::IO, t::Average)` (printing.jl L49-54)
  - Ensure Average type check works in v4
- [ ] **4.1.2** Update QTerm/QMul display (printing.jl L21-44)
  - Handle Const in arguments
- [ ] **4.1.3** Update latexify `_to_expression` (latexify_recipes.jl L134-139)
  - Handle Const in arguments
  - Verify iscall/operation/arguments usage
- [ ] **4.1.4** Test printing output

**Files to modify:**
- `/src/printing.jl`
- `/src/latexify_recipes.jl`

**Tests to run:**
- Visual inspection of output
- Latexify tests if they exist

**Acceptance Criteria:**
- Printing works correctly
- Latexification works
- No display errors

---

### Phase 5: Commutator & Specialized Operations (Week 5-6)

**Goal:** Update specialized operations like commutators.

#### 5.1 Update `commutator.jl`

**File:** `src/commutator.jl`

**Changes:**

- [ ] **5.1.1** Audit direct `.arguments` field access (L40, 49, 88, 106, 120)
  - Verify still valid in v4 or use `arguments()` function
  - Handle Const wrapping if needed
- [ ] **5.1.2** Test commutator operations

**Files to modify:**
- `/src/commutator.jl`

**Tests to run:**
- Commutator tests (if exist)
- `test/test_fock.jl`, `test/test_spin.jl` (use commutators)

**Acceptance Criteria:**
- Commutators compute correctly
- All operator tests pass

---

### Phase 6: Integration & Testing (Week 6-7)

**Goal:** Comprehensive testing and final integration.

#### 6.1 Full Test Suite

- [ ] **6.1.1** Run full test suite: `test/runtests.jl`
- [ ] **6.1.2** Fix any remaining failures
- [ ] **6.1.3** Verify code quality tests pass: `test/test_code_quality.jl`
- [ ] **6.1.4** Run benchmarks: `benchmarks/runbenchmarks.jl`
  - Compare performance vs. v3/v6 baseline
  - Investigate any regressions > 20%
- [ ] **6.1.5** Test numeric conversion: `test/test_numeric_conversion.jl`
- [ ] **6.1.6** Test cluster expansions: `test/test_cluster.jl`

**Tests:** All 13 test files

**Acceptance Criteria:**
- All tests pass
- Performance within acceptable range
- No regressions in functionality

---

#### 6.2 Documentation Updates

- [ ] **6.2.1** Update docstrings mentioning Symbolics/SymbolicUtils APIs
- [ ] **6.2.2** Update README if needed
- [ ] **6.2.3** Update migration guide in docs (if exists)
- [ ] **6.2.4** Test documentation builds: `docs/make.jl`

**Files:**
- README.md
- docs/*

**Acceptance Criteria:**
- Documentation builds
- Examples work with v7/v4

---

#### 6.3 Update Project.toml Compat

- [ ] **6.3.1** Update compat entries in `Project.toml`
  ```toml
  [compat]
  SymbolicUtils = "4"
  Symbolics = "7"
  ```
- [ ] **6.3.2** Verify all dependencies compatible with new versions
- [ ] **6.3.3** Test Project instantiation from scratch

**Files:**
- `/Project.toml`

**Acceptance Criteria:**
- Project instantiates correctly
- No dependency conflicts
- CI passes (if available)

---

### Phase 7: Release Preparation (Week 7-8)

#### 7.1 Create Migration Guide

- [ ] **7.1.1** Write MIGRATION.md documenting breaking changes for SQA users
- [ ] **7.1.2** Provide examples of updated usage patterns
- [ ] **7.1.3** List deprecated APIs (if any)

#### 7.2 Prepare Release

- [ ] **7.2.1** Update CHANGELOG.md
- [ ] **7.2.2** Update version to v0.5.0 (major version bump for breaking changes)
- [ ] **7.2.3** Tag release
- [ ] **7.2.4** Register with Julia registry

---

## V. Testing Strategy

### A. Baseline Testing (Before Migration)

1. **Full test suite on v3/v6:**
   - Record all test results
   - Identify pre-existing failures
   - Benchmark performance

2. **Document critical behaviors:**
   - Operator algebra properties
   - Average computation
   - Index summation
   - Numeric conversion

### B. Incremental Testing (During Migration)

1. **Test each phase immediately:**
   - Run relevant test files after each file modification
   - Do not proceed to next phase until current phase tests pass

2. **Regression testing:**
   - After each change, run broader test suite
   - Catch unexpected interactions early

### C. Final Validation (After Migration)

1. **Full test suite:**
   - All 13 test files must pass
   - No new failures vs. baseline

2. **Performance benchmarking:**
   - Run benchmarks/runbenchmarks.jl
   - Compare to v3/v6 baseline
   - Investigate > 20% regressions

3. **Manual testing:**
   - Try key workflows from README examples
   - Test edge cases
   - Verify printing/display

---

## VI. Risk Assessment

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| AddMul API incompatible with current pattern | Medium | High | Early prototype in Phase 0; find alternative construction |
| BasicSymbolic{T} type changes break too much | Medium | Critical | Consider wrapper types; investigate runtime type checking |
| Performance regression > 50% | Low | High | Benchmark early; work with SymbolicUtils devs if needed |
| Const wrapping breaks tree traversal | Medium | High | Comprehensive audit of arguments() calls; add unwrap_const |
| Shape field causes unexpected issues | Low | Medium | Default to scalar shape; test incrementally |
| Metadata incompatibilities | Low | Medium | Test setmetadata early; check v4 docs |
| FnType no longer supported | Medium | High | Prototype sym_average early; find alternative if broken |
| Hidden dependencies on removed APIs | Low | High | Grep codebase for all SymbolicUtils/Symbolics imports |

---

## VII. Rollback Plan

If migration proves infeasible:

1. **Partial migration:**
   - Identify minimum viable v4 compatibility
   - Keep some types on v3 patterns if possible

2. **Fork SymbolicUtils:**
   - Maintain v3 fork for SQA
   - Not recommended, high maintenance burden

3. **Delay migration:**
   - Stay on v3/v6 until ecosystem settles
   - Monitor for v4/v7 adoption issues

4. **Redesign:**
   - Consider alternative symbolic backends
   - Extreme option, very high effort

---

## VIII. Success Criteria

Migration is complete when:

- [ ] All 13 test files pass
- [ ] Performance within 20% of v3/v6 baseline
- [ ] Project.toml updated with compat entries
- [ ] Documentation updated and builds
- [ ] CHANGELOG.md updated
- [ ] Migration guide written
- [ ] Code quality checks pass (Aqua, JET, ExplicitImports)
- [ ] No deprecation warnings
- [ ] Release tagged as v0.5.0

---

## IX. Timeline Estimate

- **Phase 0:** 1 week (Preparation)
- **Phase 1:** 2 weeks (Foundation)
- **Phase 2:** 2 weeks (Indexed types)
- **Phase 3:** 1 week (Utilities)
- **Phase 4:** 1 week (Display)
- **Phase 5:** 1 week (Commutators)
- **Phase 6:** 2 weeks (Testing)
- **Phase 7:** 1 week (Release)

**Total:** ~11 weeks (~2.5 months)

**Note:** Timeline assumes single developer, part-time work. Can be accelerated with multiple contributors or full-time effort.

---

## X. Relevant File & Line References

### Critical Files (Require Immediate Attention)

| File | Lines | Critical APIs | Issue Severity |
|------|-------|---------------|----------------|
| `src/average.jl` | 11, 13-16, 19, 22, 26, 43-47 | BasicSymbolic{T}, FnType, symtype, Add/Mul | 🔴 CRITICAL |
| `src/qnumber.jl` | 51-59, 109, 114-120, 229, 240 | promote_symtype, symtype, maketerm, hash | 🔴 CRITICAL |
| `src/cnumber.jl` | 16-19, 119-123, 157-174 | Sym{T}, BasicSymbolic{T} | 🔴 CRITICAL |
| `src/index_average.jl` | 6, 40-45, 87-99, 160-217, 268-372 | BasicSymbolic{T}, Add/Mul, metadata | 🔴 CRITICAL |
| `src/indexing.jl` | 4-6, 49-86, 240-346 | Sym{T}, Ranges, tree traversal | 🟡 HIGH |

### Medium Priority Files

| File | Lines | APIs to Update | Issue Severity |
|------|-------|----------------|----------------|
| `src/index_utils.jl` | 8-28, 84-100 | BasicSymbolic, tree reconstruction | 🟡 MEDIUM |
| `src/utils.jl` | 89-96, 162-180 | Tree traversal, maketerm | 🟡 MEDIUM |
| `src/index_double_sums.jl` | 16-128 | arguments, iscall | 🟡 MEDIUM |
| `src/index_numbered_operator.jl` | 15-31 | Tree traversal | 🟡 MEDIUM |

### Low Priority Files (Minor Updates)

| File | Lines | APIs to Update | Issue Severity |
|------|-------|----------------|----------------|
| `src/printing.jl` | 21-54 | show_term, arguments | 🟢 LOW |
| `src/latexify_recipes.jl` | 134-139 | Tree traversal | 🟢 LOW |
| `src/commutator.jl` | 40, 49, 88, 106, 120 | .arguments access | 🟢 LOW |

---

## XI. Additional Resources & Links

- **SymbolicUtils v4 Docs:** https://docs.sciml.ai/SymbolicUtils/stable/
- **Symbolics v7 Docs:** https://docs.sciml.ai/Symbolics/stable/
- **Moshi.jl:** https://github.com/JuliaComputing/Moshi.jl
- **EnumX.jl (for AddMulVariant):** https://github.com/fredrikekre/EnumX.jl
- **TermInterface.jl:** https://github.com/JuliaSymbolics/TermInterface.jl
- **SciML Discourse:** https://discourse.julialang.org/c/domain/sciml/60 (for questions)

---

## XII. Questions for SymbolicUtils/Symbolics Maintainers

Before starting migration, clarify with upstream:

1. Is `FnType{Tuple{...}, RetType}` still supported in v4? (Used in `sym_average`)
2. What's the recommended way to construct AddMul variants? (AddMul(..., variant=AddMulVariant.Add)?)
3. Can `Sym{CustomType}` still be used where CustomType is user-defined? Or must CustomType be a vartype?
4. Does `setmetadata` interact with the shape field? Any gotchas?
5. Is there a migration guide for library authors (not just users)?
6. Any plans to support custom Symbolic subtypes again, or is BasicSymbolic the only path?

---

## Contact

For questions or help with this migration, tag:
- @maintainer1 (replace with actual maintainer)
- @maintainer2

Or open a discussion in: https://github.com/qojulia/SecondQuantizedAlgebra.jl/discussions
