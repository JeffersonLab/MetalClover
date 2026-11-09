#!/usr/bin/env python3
"""Print a compact read-only snapshot of the long 32^4 study."""
import csv
import io
import json
import sys
from collections import defaultdict
from pathlib import Path

p=Path(sys.argv[1] if len(sys.argv)>1 else 'results/threelevel32-compare.csv')
rows=[]
if p.exists():
    complete=''.join(s for s in p.read_text().splitlines(keepends=True) if s.endswith('\n'))
    rows=list(csv.DictReader(io.StringIO(complete)))
groups=defaultdict(list)
for r in rows:groups[r['backend']+'/'+r['method']].append(r)
tile=int(rows[0]['tile_rhs']) if rows else 4
cases={}
for key,rs in groups.items():
    v={'tiles':len(rs)}
    for rhs in (96,192):
        selected=[r for r in rs if int(r['trial'])==0 and int(r['offset'])<rhs]
        if sorted(int(r['offset']) for r in selected)==list(range(0,rhs,tile)):
            v[str(rhs)+'_seconds']=round(sum(float(r['total_ms']) for r in selected)/1000,3)
    cases[key]=v
log=p.with_suffix('.log')
last=log.read_text().splitlines()[-1] if log.exists() and log.stat().st_size else ''
print(json.dumps({'tiles':len(rows),'expected_tiles':6*192//tile,'tile_rhs':tile,
                  'worst_accepted_residual':max((float(r['max_true_residual']) for r in rows),default=0),
                  'cases':cases,'latest':last},separators=(',',':')))
