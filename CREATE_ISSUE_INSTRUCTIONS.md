# Instructions: Create GitHub Issue for Symbolics v7 Migration

## Automated Method (Recommended)

If you have the GitHub CLI (`gh`) installed and authenticated:

```bash
cd /path/to/SecondQuantizedAlgebra.jl
./.github/issues/create-migration-issue.sh
```

This will automatically create the issue with proper title, labels, and content.

## Manual Method

If you don't have `gh` CLI or prefer to create manually:

### Step 1: Navigate to GitHub
Go to: https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/new

### Step 2: Set Issue Title
```
Migration Plan: Upgrade to Symbolics v7 & SymbolicUtils v4
```

### Step 3: Copy Issue Body
Open the file `.github/issues/symbolics-v7-migration-plan.md` and copy the **entire contents** into the issue description field.

### Step 4: Add Labels
Add these labels to the issue:
- `enhancement`
- `breaking`

### Step 5: Submit
Click "Submit new issue"

## What Happens Next?

Once the issue is created, the migration plan will be available for:
1. **Team discussion** - Review and refine the plan
2. **Task assignment** - Assign specific phases to contributors  
3. **Progress tracking** - Check off tasks as they're completed
4. **Reference** - Link to specific sections when implementing changes

## Migration Plan Highlights

The issue contains:
- **986 lines** of detailed migration documentation
- **50+ actionable tasks** with checkboxes
- **13 source files** mapped with line-level references
- **7 migration phases** with acceptance criteria
- **Risk assessment** and rollback strategies
- **Testing approach** and success criteria
- **Timeline estimate**: ~11 weeks
- **Updated compat targets**: `SymbolicUtils = "4"`, `Symbolics = "7"`

## Questions?

See the `MIGRATION_PLAN_SUMMARY.md` file for a high-level overview, or open a discussion in the repository.
