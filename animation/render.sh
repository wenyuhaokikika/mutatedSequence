#!/usr/bin/env bash
# Render the explainer with Remotion inside a Linux container.
# Usage: ./render.sh [output-path]
set -euo pipefail
cd "$(dirname "$0")"
OUT=${1:-out/explainer.mp4}
mkdir -p out

if [ -n "${MSYS:-}" ] || [ "${OS:-}" = "Windows_NT" ]; then
  export MSYS_NO_PATHCONV=1
  HOST_PWD="$(pwd -W)"
else
  HOST_PWD="$PWD"
fi

docker build -f Dockerfile.render -t mutatedsequence-remotion .
docker run --rm \
  -v "$HOST_PWD":/app \
  -w /app \
  mutatedsequence-remotion \
  bash -lc "if [ ! -d node_modules/remotion ]; then npm ci; fi && \
    npx remotion render src/index.ts Explainer '$OUT' --concurrency=4"

echo "rendered: $OUT"
