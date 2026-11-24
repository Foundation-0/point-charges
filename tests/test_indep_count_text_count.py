"""A count stored as text in the input file must be compared as a number.

Found in a review: with "count": "4" the checker printed a mismatch between 4 and 4. It now agrees, and a count
field that is not a whole number is reported as a mismatch.
Usage: python3 tests/test_indep_count_text_count.py   (a few seconds; needs python-flint and NumPy)"""
import json, os, subprocess, sys, tempfile
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cert = json.load(open(os.path.join(ROOT, "certificates", "n3-triangle.json")))
ok = True
def run(count):
    c = dict(cert); c["count"] = count
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f: json.dump(c, f)
    try: return subprocess.run([sys.executable, "-I", "-O", os.path.join(ROOT, "indep_count.py"), f.name], capture_output=True, text=True)
    finally: os.unlink(f.name)
for count, want_exit, want_text in (("4", 0, "AGREES with claimed count 4"), (4, 0, "AGREES with claimed count 4"), ("5", 2, "MISMATCH: file claims 5, independent count is 4"), ("four", 2, "not a whole number")):
    r = run(count); good = r.returncode == want_exit and want_text in r.stdout; ok &= good
    print("%-8r exit %d %s" % (count, r.returncode, "OK" if good else "WRONG"))
print("ALL CASES BEHAVED" if ok else "FAILURE"); sys.exit(0 if ok else 1)
