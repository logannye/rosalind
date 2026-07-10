#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"
BIN="${ROSALIND_BIN:-target/debug/rosalind}"
if [ ! -x "$BIN" ]; then cargo build --bin rosalind; fi
D="$(mktemp -d "${TMPDIR:-/tmp}/rosalind-giab-smoke.XXXXXX")"
trap 'rm -rf "$D"' EXIT
cat > "$D/ref.fa" <<'EOF'
>chr20
AAAAAAAAAAAAAAAAAAAA
EOF
cat > "$D/truth.vcf" <<'EOF'
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr20	5	.	A	C	50	PASS	.	GT	0/1
EOF
cat > "$D/calls.vcf" <<'EOF'
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr20	5	.	A	C	50	PASS	.	GT	0/1
chr20	8	.	A	G	5	LowQual	.	GT	0/1
EOF
printf 'chr20\t0\t20\n' > "$D/confident.bed"
"$BIN" eval-germline --reference "$D/ref.fa" --calls "$D/calls.vcf" \
  --truth "$D/truth.vcf" --regions "$D/confident.bed" --calls-filter all --json > "$D/all.json"
"$BIN" eval-germline --reference "$D/ref.fa" --calls "$D/calls.vcf" \
  --truth "$D/truth.vcf" --regions "$D/confident.bed" --calls-filter pass --json > "$D/pass.json"
python3 - "$D/all.json" "$D/pass.json" <<'PY'
import json, sys
all_calls, pass_calls = (json.load(open(path)) for path in sys.argv[1:])
assert all_calls["calls_total"] == 2
assert pass_calls["calls_total"] == 1
assert pass_calls["f1"] == 1.0
PY
echo "synthetic chr20-shaped GIAB smoke: PASS"
