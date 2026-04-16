#!/usr/bin/env bash
# Local CI simulation — mirrors .github/workflows/ci.yaml
# Run before every push: bash ci-local.sh
# Or install as pre-push hook: bash ci-local.sh is called automatically.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

ok()   { echo -e "${GREEN}[OK]${NC} $*"; }
fail() { echo -e "${RED}[FAIL]${NC} $*"; exit 1; }
step() { echo -e "\n${YELLOW}==>${NC} $*"; }

cd "$(git rev-parse --show-toplevel)"

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

# ── 2. testthat ───────────────────────────────────────────────────────────────
step "Running tests (testthat)"
Rscript -e '
  # Source all R files in dependency order
  files <- c("R/alphas.R", "R/gram.R", "R/samplers.R", "R/domain.R",
             "R/utils.R", "R/core.R", "R/density_lp.R", "R/methods.R")
  invisible(lapply(files, source))
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
