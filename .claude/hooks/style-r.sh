#!/usr/bin/env bash
# PostToolUse hook (Edit|Write): auto-style touched .R files with styler.
# Receives the tool-call JSON on stdin; always exits 0 so styling problems
# never block the agent — lintr in CI is the enforcing gate.

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
  # generated files are never styled (and must never be hand-edited)
  */R/RcppExports.R | */R/stanmodels.R | */src/*) exit 0 ;;
  *.R | *.r)
    Rscript --no-init-file -e \
      "if (requireNamespace('styler', quietly = TRUE)) styler::style_file('$file')" \
      >/dev/null 2>&1
    ;;
esac

exit 0
