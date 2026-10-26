#!/usr/bin/env python3
"""Resolve approved PbPb B preprocessed DATA. No raw-input fallback or plotting."""
import argparse
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = REPO.parent / 'XGBoost/input/B_pbpb_precut_data.manifest.json'
SELECTION = '(Bpt > 7.5) and (abs(By) < 2.4) and (Bchi2Prob > 0.05) and (Bnorm_svpvDistance_2D > 2) and (Bmass > 5) and (Bmass < 6)'
ROOT_SELECTION = '(Bpt > 7.5) && (abs(By) < 2.4) && (Bchi2Prob > 0.05) && (Bnorm_svpvDistance_2D > 2) && (Bmass > 5) && (Bmass < 6)'
TREES = {'Bu': 'ntKp', 'Bs': 'ntphi', 'Bd': 'ntKstar'}
SIDEBANDS = {'Bu': ((5.03, 5.13), (5.43, 5.53)), 'Bd': ((5.03, 5.13), (5.43, 5.53)), 'Bs': ((5.02, 5.22), (5.52, 5.72))}
REQUIRED = ('Bmass', 'Bpt', 'By', 'Bnorm_svpvDistance_2D', 'Bchi2Prob')

def in_sideband(channel, mass):
    """Strict open intervals; supports scalar or NumPy array input."""
    (a, b), (c, d) = SIDEBANDS[channel]
    return ((mass > a) & (mass < b)) | ((mass > c) & (mass < d))

def resolve_input(channel, year, *, scope, required_branches=(), manifest=DEFAULT_MANIFEST):
    """Validate metadata, returning a provenance-ready plotting input contract.

    scope='full' retains the full preprocessed mass window; 'sideband' adds
    the approved sideband OR. Neither scope can recover removed candidates.
    Missing or stale files raise errors; there is no fallback to raw DATA.
    """
    import uproot
    if channel not in TREES or year not in (2023, 2024):
        raise ValueError('Only Bu/Bs/Bd PbPb23/24 DATA are supported')
    if scope not in ('full', 'sideband'):
        raise ValueError('Explicit scope must be full or sideband')
    manifest = Path(manifest)
    payload = manifest.read_bytes()
    m = json.loads(payload)
    if m['status'] != 'PASS' or m['selection'] != SELECTION:
        raise ValueError('Unaccepted or changed preprocessing contract; review required')
    e = m['files'][f'{channel}_{year}']
    p = Path(e['path'])
    if e['tree'] != TREES[channel] or e['output']['path'] != str(p):
        raise ValueError('Manifest path/tree mismatch')
    if len(e['sha256']) != 64 or any(c not in '0123456789abcdef' for c in e['sha256']):
        raise ValueError('Invalid producer SHA256')
    before = p.stat()
    if (before.st_size, before.st_mtime_ns) != (e['output']['size_bytes'], e['output']['mtime_ns']):
        raise ValueError('Input fingerprint changed')
    with uproot.open(p) as f:
        t = f[e['tree']]
        if t.num_entries != e['entries']:
            raise ValueError('Entry count mismatch')
        missing = set(REQUIRED).union(required_branches) - set(t.keys())
        if missing:
            raise ValueError(f'Missing branches: {sorted(missing)}')
        branch_types = t.typenames()
    after = p.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('Input changed during validation')
    intervals = SIDEBANDS[channel]
    if not all(5 < lo < hi < 6 for lo, hi in intervals):
        raise ValueError('Sideband extends outside preprocessed mass range')
    sideband = ' || '.join(f'((Bmass > {lo:.2f}) && (Bmass < {hi:.2f}))' for lo, hi in intervals)
    return {
        'channel': channel, 'year': year, 'kind': 'DATA', 'path': str(p),
        'tree': e['tree'], 'entries': e['entries'], 'output': e['output'],
        'producer_sha256': e['sha256'], 'source_sha256_recomputed': False,
        'manifest': str(manifest), 'manifest_sha256': hashlib.sha256(payload).hexdigest(),
        'preprocessing_selection': SELECTION, 'preprocessing_selection_ROOT': ROOT_SELECTION,
        'scope': scope, 'additional_selection_ROOT': sideband if scope == 'sideband' else '1',
        'sideband_intervals_GeV': intervals, 'sideband_selection_ROOT': sideband,
        'sideband_applied': scope == 'sideband', 'sideband_definition': 'ML request 2026-10-14; strict open intervals',
        'branch_types': branch_types, 'validation': 'metadata PASS; producer SHA recorded, not recomputed',
    }

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--channel', choices=list(TREES), required=True)
    parser.add_argument('--year', type=int, choices=[2023, 2024], required=True)
    parser.add_argument('--scope', choices=['full', 'sideband'], required=True)
    parser.add_argument('--require-branch', action='append', default=[])
    args = parser.parse_args()
    print(json.dumps(resolve_input(args.channel, args.year, scope=args.scope,
                                  required_branches=args.require_branch), indent=2))

if __name__ == '__main__':
    main()
