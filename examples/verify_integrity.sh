#!/usr/bin/env sh
# Check every file of the repository against its recorded SHA-256.
set -e
sha256sum -c manifests/SHA256SUMS --quiet && echo "all files match manifests/SHA256SUMS"
