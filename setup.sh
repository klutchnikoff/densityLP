#!/usr/bin/env bash
# Run once after cloning: bash setup.sh
# Configures git to use the versioned hooks in .githooks/

set -euo pipefail

cd "$(git rev-parse --show-toplevel)"
git config core.hooksPath .githooks
echo "Git hooks configured (.githooks/)."
echo "  pre-commit : docs + format"
echo "  pre-push   : docs + format + tests + R CMD check"
