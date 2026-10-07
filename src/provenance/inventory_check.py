"""Static evidence lint, not scientific-tool certification or API execution."""
import json
import math
import sys
from pathlib import Path

def nonfinite_paths(value, path='$'):
    if isinstance(value, float) and not math.isfinite(value):
        return [path]
    if isinstance(value, dict):
        return [p for key, item in value.items() for p in nonfinite_paths(item, path + '.' + key)]
    if isinstance(value, list):
        return [p for i, item in enumerate(value) for p in nonfinite_paths(item, path + '[' + str(i) + ']')]
    return []

def inspect_inventory(doc):
    issues = []
    entries = doc.get('tools')
    if not isinstance(entries, list):
        return {'certified': False, 'issues': ['tools must be a list'], 'recorded_entries': 0}
    names = [t.get('tool') for t in entries if isinstance(t, dict)]
    if len(names) != len(entries) or any(not isinstance(n, str) or not n for n in names):
        issues.append('malformed tool entry or name')
    valid_names = [n for n in names if isinstance(n, str)]
    if len(set(valid_names)) != len(valid_names):
        issues.append('duplicate tool names')
    actual_ok = sum(isinstance(t, dict) and t.get('status') == 'ok' for t in entries)
    if doc.get('n_tools_ok') != actual_ok:
        issues.append('reported ok count does not match entries')
    if doc.get('n_tools_attempted') != len(entries):
        issues.append('reported attempted count does not match entries')
    for i, t in enumerate(entries):
        if not isinstance(t, dict):
            continue
        label = t.get('tool', str(i))
        if t.get('status') not in ('ok', 'failed'):
            issues.append(str(label) + ': unknown status')
        if t.get('status') == 'ok' and not isinstance(t.get('result'), dict):
            issues.append(str(label) + ': ok without result object')
        for p in nonfinite_paths(t):
            issues.append(str(label) + ': non-finite value at ' + p)
        if t.get('status') == 'failed' and not t.get('error'):
            issues.append(str(label) + ': failed without error')
    for k in ('code_revision', 'data_sha256', 'package_versions', 'executed_at'):
        if not doc.get(k):
            issues.append('missing execution provenance: ' + k)
    return {'certified': False, 'recorded_entries': len(entries), 'recorded_ok': actual_ok,
            'issues': issues, 'boundary': 'Passing structural checks cannot certify source identity, valid analysis, distinct tools, novelty or execution.'}

if __name__ == '__main__':
    p = Path(sys.argv[1])
    # Historical Python JSON may contain Infinity; keep it visible for diagnosis.
    print(json.dumps(inspect_inventory(json.loads(p.read_text())), indent=2, allow_nan=False))
