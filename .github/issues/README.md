# GitHub Issues - Migration Plans

This directory contains detailed migration plans that should be opened as GitHub issues.

## Symbolics v7 & SymbolicUtils v4 Migration

**File:** `symbolics-v7-migration-plan.md`

**To create the GitHub issue:**

1. Navigate to: https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/new
2. Set title: `Migration Plan: Upgrade to Symbolics v7 & SymbolicUtils v4`
3. Copy the entire contents of `symbolics-v7-migration-plan.md` into the issue body
4. Add labels: `enhancement`, `breaking`
5. Submit the issue

Alternatively, if you have the GitHub CLI (`gh`) installed with authentication:

```bash
cd /path/to/SecondQuantizedAlgebra.jl
gh issue create --title "Migration Plan: Upgrade to Symbolics v7 & SymbolicUtils v4" \
  --body-file .github/issues/symbolics-v7-migration-plan.md \
  --label "enhancement,breaking"
```
