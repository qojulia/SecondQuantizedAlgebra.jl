#!/bin/bash
# Script to create the Symbolics v7 migration GitHub issue
# This requires the GitHub CLI (gh) to be installed and authenticated

set -e

REPO="qojulia/SecondQuantizedAlgebra.jl"
TITLE="Migration Plan: Upgrade to Symbolics v7 & SymbolicUtils v4"
ISSUE_FILE=".github/issues/symbolics-v7-migration-plan.md"

echo "Creating GitHub issue: $TITLE"
echo "Repository: $REPO"
echo "Issue body from: $ISSUE_FILE"
echo ""

if ! command -v gh &> /dev/null; then
    echo "ERROR: GitHub CLI (gh) not found. Please install it from https://cli.github.com/"
    exit 1
fi

if [ ! -f "$ISSUE_FILE" ]; then
    echo "ERROR: Issue file not found: $ISSUE_FILE"
    exit 1
fi

echo "Creating issue..."
gh issue create \
    --repo "$REPO" \
    --title "$TITLE" \
    --body-file "$ISSUE_FILE" \
    --label "enhancement,breaking"

echo ""
echo "Issue created successfully!"
