#!/usr/bin/env python3
"""Check an example run against the committed hashes in docs/expected-example.json.

    python scripts/regression_compare.py --check-expected docs/expected-example.json --new-out DIR

DIR is the output_dir of a run of the example manifest. Each arm packet (sources/NNNN/<arm>.json) is hashed
with the timing field removed and compared with the committed digest. Exit status 0 if all match.
"""
import argparse, hashlib, json, pathlib, sys


def packet_digest(path):
    """sha256 of an arm packet with the timing field removed (canonical JSON), comparable across machines."""
    packet = {k: v for k, v in json.loads(pathlib.Path(path).read_text()).items() if k != 'seconds'}
    return hashlib.sha256(json.dumps(packet, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--check-expected', required=True, help='committed expected-example.json')
    p.add_argument('--new-out', required=True, help='output_dir of an example run')
    args = p.parse_args()
    exp = json.loads(pathlib.Path(args.check_expected).read_text())
    bad = []
    for rank, arms in exp['arm_packet_sha256'].items():
        for arm, want in arms.items():
            got = packet_digest(pathlib.Path(args.new_out) / 'sources' / f'{int(rank):04d}' / f'{arm}.json')
            if got != want:
                bad.append(f'rank {rank} {arm}: {got} != {want}')
    print(json.dumps(dict(checked=sum(len(a) for a in exp['arm_packet_sha256'].values()), mismatches=bad, passed=not bad), indent=2))
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main())
