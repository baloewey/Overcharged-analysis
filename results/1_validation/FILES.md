# 1_validation

Is the method right? Custom anchoring, its resolution limit, new wall shapes, tracker, seeds.

Copied by `tools/build_results.py`; regenerate rather than edit. Files marked KEY are the main figures; SUPERSEDED files are kept for the record only.

| File | What it is | Copied from |
|---|---|---|
| `resolution_limit.txt` | Anchoring resolution limit: clean to rho = 0.15, fails from 0.25 (passive and active alike) | `runs/cavity/resolution_scan/resolution_results.txt` |
| `resolution_limit_active_checks.txt` | Per-run audit, active resolution scan | `runs/cavity/resolution_scan/act_check.txt` |
| `resolution_limit_passive_checks.txt` | Per-run audit, passive resolution scan | `runs/cavity/resolution_scan/pas_check.txt` |
| `cassini_walls_passive_tests.png` | Cassini wall implementation: a = 0 reproduces the circle; passive ovals with a concave neck | `runs/cassini/tests_T1T2.png` |
| `cassini_passive_pinched_oval.png` | Passive pinched ovals keep all 2k defects at the wall (13-cell neck) | `runs/cassini/P2/passive_check.png` |
| `server_binary_validation_cassini.png` | Cassini binary on the server = office build (2026-10-03) | `runs/cassini/server_validation/check.png` |
| `server_benchmark.txt` | Server throughput benchmark | `runs/server_benchmark.txt` |
