#!/bin/bash
# prepare_docs.sh — prepare example notebooks and build Sphinx HTML docs
#
# Usage:
#   bash bin/prepare_docs.sh          # prepare examples only
#   bash bin/prepare_docs.sh build    # prepare examples and build HTML

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$(dirname "$DIR")"

echo "Preparing example notebooks for documentation..."
python3 "$DIR/examples_doc.py"

if [[ "${1:-}" == "build" ]]; then
    echo "Building Sphinx HTML documentation..."
    cd "$ROOT_DIR/docs"
    make html
    echo "Documentation built at docs/_build/html/index.html"
fi

echo "Done."
