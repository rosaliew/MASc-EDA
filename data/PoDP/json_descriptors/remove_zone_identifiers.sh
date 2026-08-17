#!/usr/bin/env bash
# Remove Windows Zone.Identifier alternate-data-stream files from this folder.
set -euo pipefail
find "$(dirname "$0")" -iname "*Zone.Identifier*" -type f -print -delete
