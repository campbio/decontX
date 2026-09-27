#!/usr/bin/env bash
# SessionStart hook for Claude Code: prints the development standards so
# Claude Code adds them to the session's context.
#
# Copy this file to dev/hooks/load-standards.sh in each package and register
# it in .claude/settings.json (see ADOPTING.md).
#
# Source: https://github.com/campbio/r-bioc-dev-standards
# If you use a fork, change the default URL below to point at it.
#
# To test an edit on a branch for one session:
#   R_BIOC_STANDARDS_URL=https://raw.githubusercontent.com/<owner>/r-bioc-dev-standards/<branch>/standards.md claude
#
# If GitHub is unreachable, the last downloaded copy is used instead.
# Always exits 0 so a network problem never blocks a session.

set -u

url="${R_BIOC_STANDARDS_URL:-https://raw.githubusercontent.com/campbio/r-bioc-dev-standards/main/standards.md}"
cache_dir="${XDG_CACHE_HOME:-$HOME/.cache}/r-bioc-dev-standards"
cache="$cache_dir/standards.md"

mkdir -p "$cache_dir" 2>/dev/null
tmp="$(mktemp "$cache_dir/download.XXXXXX" 2>/dev/null)" || tmp=""

fresh=0
if [ -n "$tmp" ] && curl -fsSL --max-time 10 "$url" -o "$tmp" 2>/dev/null \
   && [ -s "$tmp" ]; then
  mv -f "$tmp" "$cache" && fresh=1
fi
[ -n "$tmp" ] && [ -f "$tmp" ] && rm -f "$tmp"

if [ -s "$cache" ]; then
  if [ "$fresh" -eq 0 ]; then
    echo "NOTE: The development standards could not be downloaded, so a"
    echo "cached copy is being used and may be out of date. Tell the user."
    echo
  fi
  cat "$cache"
else
  echo "WARNING: The development standards (r-bioc-dev-standards) could not"
  echo "be loaded and no cached copy exists. Tell the user before starting"
  echo "any work."
fi
exit 0
