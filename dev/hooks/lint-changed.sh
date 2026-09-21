#!/usr/bin/env bash
# PostToolUse hook (Edit|Write): lint a touched .R file and report
# findings back to the agent. REPORT-ONLY — it never rewrites files:
# with a large pre-existing style backlog, auto-styling every touched
# file would bury real changes in formatting diffs. Styling is applied
# deliberately, per file, by a human-approved change.
#
# Exit 2 feeds stderr back to the agent as feedback; the edit itself is
# never blocked (the tool has already run).

set -u

input=$(cat)

if command -v jq >/dev/null 2>&1; then
  file=$(printf '%s' "$input" | jq -r '.tool_input.file_path // empty')
else
  file=$(printf '%s' "$input" |
    sed -n 's/.*"file_path"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' |
    head -n 1)
fi

[ -n "${file:-}" ] || exit 0
[ -f "$file" ] || exit 0

case "$file" in
  # generated files are never linted (and must never be hand-edited)
  */R/RcppExports.R | */R/stanmodels.R | */src/*) exit 0 ;;
  *.R | *.r) ;;
  *) exit 0 ;;
esac

# pass the path via the environment, not string interpolation, so paths
# with quotes cannot break out of (or into) the R expression
lints=$(LINT_FILE="$file" Rscript --no-init-file -e \
  "if (requireNamespace('lintr', quietly = TRUE)) { l <- lintr::lint(Sys.getenv('LINT_FILE')); if (length(l) > 0) print(l) }" \
  2>/dev/null)

if [ -n "$lints" ]; then
  {
    printf '%s\n' "$lints" | head -n 40
    echo "lintr findings in $file (first 40 lines shown)."
    echo "Fix any your change introduced; leave pre-existing debt alone."
  } >&2
  exit 2
fi

exit 0
