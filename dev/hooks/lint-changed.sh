#!/usr/bin/env bash
# PostToolUse hook for Claude Code: lints an R file right after Claude edits
# it and passes the results back to Claude. It reports only and never
# rewrites the file: auto-styling on every edit would bury the real change in
# formatting noise.
#
# Copy this file to dev/hooks/lint-changed.sh in each package and register
# it in .claude/settings.json (see ADOPTING.md).
#
# Source: https://github.com/campbio/r-bioc-dev-standards
#
# Claude Code sends the hook payload as JSON on stdin. Plain stdout from a
# PostToolUse hook only reaches the debug log, so the lints are returned as
# JSON in hookSpecificOutput.additionalContext, which Claude does see.
# Always exits 0: a lint is information, not a reason to fail the edit.

set -uo pipefail

payload="$(cat)"

# Extract .tool_input.file_path without depending on jq.
if command -v jq > /dev/null 2>&1; then
  file="$(printf '%s' "$payload" \
    | jq -r '.tool_input.file_path // empty' 2> /dev/null)"
elif command -v python3 > /dev/null 2>&1; then
  file="$(printf '%s' "$payload" | python3 -c '
import json, sys
try:
    print(json.load(sys.stdin).get("tool_input", {}).get("file_path", ""))
except Exception:
    print("")')"
else
  exit 0
fi

[ -n "$file" ] && [ -f "$file" ] || exit 0
case "$file" in
  *.R|*.r) ;;
  *) exit 0 ;;
esac
command -v Rscript > /dev/null 2>&1 || exit 0

# One line per lint, at most 25, so a file with a lint backlog can't flood
# the session.
lints="$(Rscript -e '
  f <- commandArgs(trailingOnly = TRUE)[1]
  if (!requireNamespace("lintr", quietly = TRUE)) quit(status = 0)
  l <- tryCatch(as.data.frame(lintr::lint(f)), error = function(e) NULL)
  if (is.null(l) || nrow(l) == 0) quit(status = 0)
  n <- nrow(l)
  l <- utils::head(l, 25)
  cat(sprintf("lintr found %d lint(s) in %s. Fix those on lines you changed;",
              n, basename(f)),
      "leave existing lints elsewhere in the file for a separate PR.\n")
  cat(sprintf("  line %d: %s [%s]\n", l$line_number, l$message, l$linter),
      sep = "")
  if (n > 25) cat(sprintf("  ... and %d more\n", n - 25))
' "$file" 2> /dev/null)"

[ -n "$lints" ] || exit 0

if command -v jq > /dev/null 2>&1; then
  jq -n --arg ctx "$lints" \
    '{hookSpecificOutput: {hookEventName: "PostToolUse", additionalContext: $ctx}}'
elif command -v python3 > /dev/null 2>&1; then
  printf '%s' "$lints" | python3 -c '
import json, sys
print(json.dumps({"hookSpecificOutput": {
    "hookEventName": "PostToolUse", "additionalContext": sys.stdin.read()}}))'
fi
exit 0
