#!/usr/bin/env bash
# Export the full transcript of the current opencode session as JSON,
# compress it with xz -9, and stage it for commit. Run before every push
# (see AGENTS.md "Session transcripts").
set -euo pipefail

repo="$(git rev-parse --show-toplevel)"
db="$HOME/.local/share/opencode/opencode.db"
out_dir="$repo/transcripts"
mkdir -p "$out_dir"

session_id="$(python3 - "$db" "$repo" <<'EOF'
import sqlite3, sys
db, repo = sys.argv[1], sys.argv[2]
con = sqlite3.connect(db)
row = con.execute(
    "SELECT id FROM session WHERE directory = ? AND parent_id IS NULL "
    "ORDER BY time_updated DESC LIMIT 1",
    (repo,),
).fetchone()
if row is None:
    sys.exit("no opencode session found for " + repo)
print(row[0])
EOF
)"

json="$out_dir/$session_id.json"
xz_file="$out_dir/$session_id.json.xz"

opencode export "$session_id" > "$json"
python3 - "$json" <<'EOF'
import json, sys
with open(sys.argv[1]) as f:
    data = json.load(f)
if not data.get("messages"):
    sys.exit("export has no messages: " + sys.argv[1])
EOF
xz -9 -f "$json"
git add "$xz_file"
echo "staged $xz_file ($(stat -c%s "$xz_file") bytes, session $session_id)"