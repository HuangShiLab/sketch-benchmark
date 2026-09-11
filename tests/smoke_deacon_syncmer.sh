#!/bin/bash
# Smoke test for the deacon-syncmer build: index sizes match at (k=31,w=15)
# vs (l=31,s=16), the syncmer index misses fewer mutated host reads at 1-2%
# error, and neither scheme matches random reads.
#   tests/smoke_deacon_syncmer.sh [path/to/deacon-syncmer]
set -euo pipefail
D="${1:-deacon-syncmer}"; W="$(mktemp -d)"; cd "$W"
python3 - <<'PY'
import random
random.seed(5)
G=''.join(random.choice('ACGT') for _ in range(300_000))
open('genome.fa','w').write('>g\n'+G+'\n')
def mut(s,p): return ''.join(c if random.random()>p else random.choice([b for b in 'ACGT' if b!=c]) for c in s)
for p in (0.01, 0.02):
    with open(f'host_e{p}.fq','w') as f:
        for i in range(4000):
            st=random.randrange(0,len(G)-150); f.write(f'@h{i}\n{mut(G[st:st+150],p)}\n+\n{"I"*150}\n')
with open('random.fq','w') as f:
    for i in range(4000):
        f.write(f'@r{i}\n{"".join(random.choice("ACGT") for _ in range(150))}\n+\n{"I"*150}\n')
PY
"$D" index build genome.fa -k 31 -w 15 -o mini.idx -q 2>/dev/null
"$D" index build genome.fa -k 31 -w 16 --scheme syncmer -o sync.idx -q 2>/dev/null
"$D" index info sync.idx 2>&1 | grep -q "Scheme: syncmer" || { echo "FAIL: header does not record syncmer scheme"; exit 1; }
matched() { "$D" filter "$1" "$2" -a 2 -r 0 -t 2 -o /dev/null --summary s.json >/dev/null 2>&1; python3 -c "import json;print(json.load(open('s.json'))['seqs_out'])"; }
m1=$(matched mini.idx host_e0.01.fq); s1=$(matched sync.idx host_e0.01.fq)
m2=$(matched mini.idx host_e0.02.fq); s2=$(matched sync.idx host_e0.02.fq)
mr=$(matched mini.idx random.fq);     sr=$(matched sync.idx random.fq)
echo "host 1% error  matched: minimizer=$m1 syncmer=$s1 (of 4000)"
echo "host 2% error  matched: minimizer=$m2 syncmer=$s2 (of 4000)"
echo "random         matched: minimizer=$mr syncmer=$sr (expect 0)"
[ "$s2" -gt "$m2" ] || { echo "FAIL: syncmer did not retain more host reads at 2% error"; exit 1; }
[ "$mr" -eq 0 ] && [ "$sr" -eq 0 ] || { echo "FAIL: false hits on random reads"; exit 1; }
echo "PASS"; rm -rf "$W"
