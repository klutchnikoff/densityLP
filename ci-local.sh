#!/usr/bin/env bash
# Local CI simulation — mirrors .github/workflows/ci.yaml
#
# Usage:
#   bash ci-local.sh          — full CI (docs + format + tests + R CMD check)
#   bash ci-local.sh --fast   — fast checks only (docs + format), for pre-commit

set -euo pipefail

FAST=false
[[ "${1:-}" == "--fast" ]] && FAST=true

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

ok()   { echo -e "${GREEN}[OK]${NC} $*"; }
fail() { echo -e "${RED}[FAIL]${NC} $*"; exit 1; }
step() { echo -e "\n${YELLOW}==>${NC} $*"; }

cd "$(git rev-parse --show-toplevel)"

# ── 0. Docs (roxygen2) ────────────────────────────────────────────────────────
step "Checking docs are up to date (roxygen2)"
Rscript -e "roxygen2::roxygenise()" 2>/dev/null
if ! git diff --exit-code --quiet -- NAMESPACE man/; then
  echo ""
  echo "The following generated files are out of date:"
  git diff --name-only -- NAMESPACE man/
  echo ""
  echo "Run roxygen2::roxygenise(), re-stage, and commit."
  fail "Docs out of date"
fi
ok "Docs"

# ── 1. Formatting (Air) ────────────────────────────────────────────────────────
step "Checking formatting (Air)"
if ! command -v air &>/dev/null; then
  fail "air not found. Install from https://github.com/posit-dev/air"
fi
if ! air format --check .; then
  echo ""
  echo "Run 'air format .' to fix formatting, then re-stage and commit."
  fail "Formatting check failed"
fi
ok "Formatting"

if $FAST; then
  echo -e "\n${GREEN}Fast checks passed — safe to commit.${NC}"
  exit 0
fi

# ── 2. testthat ───────────────────────────────────────────────────────────────
step "Running tests (testthat)"
Rscript -e '
  install.packages(".", repos = NULL, type = "source", quiet = TRUE)
  library(densityLP)
  attach(asNamespace("densityLP"))
  library(testthat)
  results <- test_dir("tests/testthat", reporter = "progress")
  if (any(as.data.frame(results)$failed > 0)) quit(status = 1)
' || fail "Tests failed"
ok "Tests"

# ── 3. R CMD check ────────────────────────────────────────────────────────────
step "R CMD check --as-cran"
PKG=$(Rscript -e 'cat(read.dcf("DESCRIPTION")[,"Package"])')
Rscript -e 'install.packages("rcmdcheck", repos = "https://cloud.r-project.org", quiet = TRUE)' 2>/dev/null
_R_CHECK_FORCE_SUGGESTS_=false Rscript -e '
  rcmdcheck::rcmdcheck(
    args    = c("--no-manual", "--as-cran"),
    error_on = "error"
  )
' || fail "R CMD check failed"
ok "R CMD check"

echo -e "\n${GREEN}All checks passed — safe to push.${NC}"
