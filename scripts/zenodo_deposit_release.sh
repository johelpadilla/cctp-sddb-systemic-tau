#!/usr/bin/env bash
# Manual Zenodo deposition of a GitHub release archive (fallback if GitHub↔Zenodo
# webhook is not yet enabled).
#
# Prerequisites:
#   export ZENODO_TOKEN=...   # Personal access token from
#                             # https://zenodo.org/account/settings/applications/
#   curl, jq, python3
#
# Usage:
#   ./scripts/zenodo_deposit_release.sh v1.0.0
#   ./scripts/zenodo_deposit_release.sh v1.0.0 --publish
set -euo pipefail

TAG="${1:-v1.0.0}"
PUBLISH=0
[[ "${2:-}" == "--publish" ]] && PUBLISH=1

if [[ -z "${ZENODO_TOKEN:-}" ]]; then
  echo "ERROR: set ZENODO_TOKEN (Zenodo personal access token with deposit:write)."
  echo "Create at: https://zenodo.org/account/settings/applications/"
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"
API="https://zenodo.org/api"
AUTH="Authorization: Bearer ${ZENODO_TOKEN}"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

ARCHIVE="${TMP}/cctp-sddb-systemic-tau-${TAG}.zip"
echo "Downloading GitHub archive for ${TAG}..."
curl -fsSL -L \
  "https://github.com/johelpadilla/cctp-sddb-systemic-tau/archive/refs/tags/${TAG}.zip" \
  -o "$ARCHIVE"

echo "Creating Zenodo deposition..."
DEPOSIT=$(curl -fsS -H "$AUTH" -H "Content-Type: application/json" \
  -X POST "${API}/deposit/depositions" -d '{}')
DEP_ID=$(echo "$DEPOSIT" | python3 -c "import sys,json; print(json.load(sys.stdin)['id'])")
BUCKET=$(echo "$DEPOSIT" | python3 -c "import sys,json; print(json.load(sys.stdin)['links']['bucket'])")
echo "Deposition id: ${DEP_ID}"

echo "Uploading archive..."
curl -fsS -H "$AUTH" -H "Content-Type: application/octet-stream" \
  --upload-file "$ARCHIVE" \
  "${BUCKET}/cctp-sddb-systemic-tau-${TAG}.zip" >/dev/null

META_FILE="${TMP}/metadata.json"
export ZENODO_VERSION="${TAG#v}"
python3 - <<PY >"$META_FILE"
import json, os
from pathlib import Path
z = json.loads(Path(".zenodo.json").read_text())
ver = os.environ.get("ZENODO_VERSION", "1.1.0")
meta = {
    "metadata": {
        "title": z["title"],
        "upload_type": z.get("upload_type", "software"),
        "description": z["description"],
        "creators": [
            {
                "name": c["name"],
                "affiliation": c.get("affiliation", ""),
                "orcid": c.get("orcid", "").replace("https://orcid.org/", ""),
            }
            for c in z["creators"]
        ],
        "keywords": z.get("keywords", []),
        "license": z.get("license", "mit"),
        "access_right": z.get("access_right", "open"),
        "language": z.get("language", "eng"),
        "related_identifiers": z.get("related_identifiers", []),
        "version": ver,
    }
}
print(json.dumps(meta))
PY

echo "Updating metadata from .zenodo.json..."
curl -fsS -H "$AUTH" -H "Content-Type: application/json" \
  -X PUT "${API}/deposit/depositions/${DEP_ID}" \
  -d @"$META_FILE" >/dev/null

echo "Draft ready: https://zenodo.org/deposit/${DEP_ID}"

if [[ "$PUBLISH" -eq 1 ]]; then
  echo "Publishing..."
  PUB=$(curl -fsS -H "$AUTH" -H "Content-Type: application/json" \
    -X POST "${API}/deposit/depositions/${DEP_ID}/actions/publish")
  DOI=$(echo "$PUB" | python3 -c "import sys,json; print(json.load(sys.stdin).get('doi',''))")
  echo "Published DOI: ${DOI}"
  echo "Record: https://doi.org/${DOI}"
else
  echo
  echo "Review the draft in the browser, then either:"
  echo "  1) Click Publish on Zenodo, or"
  echo "  2) Re-run: $0 ${TAG} --publish"
fi
