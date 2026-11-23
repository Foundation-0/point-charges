#!/usr/bin/env python3
"""
Standalone CLI for Track S21d: Point Charges n=7 (#12) & n=8 (#5) Maxwell Frontier Certifier.

Zero floats, zero bare `assert`s, strict Zone-B fail-closed validation.
Supports `--standalone --json`, `--write [PATH]`, and `--check PATH` (defaulting to full
pure-Z/Q 384-bit Krawczyk + Bareiss Sylvester-inertia replay across all 127 equilibria).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Sequence

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from point_charges_n7_n8_oracle import (
    PointChargeFrontierError,
    build_certificate,
    canonical_json_bytes,
    load_safe_bytes,
    parse_strict_json_text,
    require_safe_path,
    require_str,
    safe_write_bytes,
    validate_certificate_payload,
)

DEFAULT_CERT_PATH = _HERE / "point_charges_n7_n8_certificate.json"


def verify_certificate_file(
    cert_path: Path,
    *,
    full_replay: bool = True,
) -> dict[str, object]:
    """Fail-closed Zone-B verification of a canonical LF-only Track S21d certificate file."""
    raw = load_safe_bytes(cert_path, "certificate_file")
    if b"\r" in raw:
        raise PointChargeFrontierError("CR_FORBIDDEN", "certificate file must use LF endings only")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise PointChargeFrontierError("INVALID_UTF8", "certificate file is not valid UTF-8") from exc
    obj = parse_strict_json_text(text, require_object=True)
    expected_raw = canonical_json_bytes(obj) + b"\n"
    if raw != expected_raw:
        raise PointChargeFrontierError("NON_CANONICAL_CERT_BYTES", "certificate file bytes are not canonical JSON + LF")
    res = validate_certificate_payload(obj, full_replay=full_replay)
    res["file_sha256"] = hashlib.sha256(raw).hexdigest()
    return res


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Track S21d Point Charges n=7 & n=8 Pure-Z/Q Krawczyk Certifier CLI")
    mode = parser.add_mutually_exclusive_group(required=False)
    mode.add_argument("--standalone", action="store_true", help="Certify all 4 configs (127 equilibria) and emit canonical certificate JSON")
    mode.add_argument("--write", nargs="?", const=str(DEFAULT_CERT_PATH), metavar="OUT_PATH", help="Write canonical certificate JSON to file")
    mode.add_argument("--check", metavar="CERT_PATH", help="Verify canonical certificate JSON (full pure-Z/Q Krawczyk replay by default)")
    parser.add_argument("--json", action="store_true", help="Emit canonical JSON output")
    parser.add_argument("--no-replay", action="store_true", help="With --check, verify schema, digests, and arithmetic consistency without full Krawczyk replay")

    try:
        if argv is not None:
            if type(argv) not in (list, tuple):
                raise PointChargeFrontierError("INVALID_CLI_ARGS", "argv must be a strict list or tuple")
            raw_argv: list[str] | None = [require_str(x, "argv", min_len=0) for x in argv]
        else:
            raw_argv = None
        args = parser.parse_args(raw_argv)

        if args.no_replay and (args.standalone or args.write is not None):
            raise PointChargeFrontierError("INVALID_CLI_ARGS", "--no-replay is valid only with --check")

        if args.standalone:
            cert = build_certificate()
            sys.stdout.write(canonical_json_bytes(cert).decode("ascii") + "\n")
            return 0

        if args.write is not None:
            out_path = require_safe_path(Path(args.write), "write_path")
            cert = build_certificate()
            raw = canonical_json_bytes(cert) + b"\n"
            safe_write_bytes(out_path, raw, "write_path")
            summary = {
                "written": True,
                "file_sha256": hashlib.sha256(raw).hexdigest(),
                "payload_sha256": cert["payload_sha256"],
                "certificate_sha256": cert["certificate_sha256"],
                "total_equilibria_certified": cert["payload"]["frontier_verdicts"]["total_equilibria_certified"],
                "issue5_n8_51_maxwell_counterexample_certified": cert["payload"]["frontier_verdicts"]["issue5_n8_51_maxwell_counterexample_certified"],
            }
            print(json.dumps(summary, sort_keys=True, separators=(",", ":")))
            return 0

        target_cert = args.check if args.check is not None else str(DEFAULT_CERT_PATH)
        cert_path = require_safe_path(Path(target_cert), "check_path")
        res = verify_certificate_file(cert_path, full_replay=not args.no_replay)
        print(json.dumps(res, sort_keys=True, separators=(",", ":")))
        return 0
    except PointChargeFrontierError as exc:
        print(json.dumps({"valid": False, "error_code": exc.code, "error": exc.message}, sort_keys=True), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
