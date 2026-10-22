#!/usr/bin/env python3
"""Turn a passive checkpoint into an active restart: switch.py IN OUT TOTAL_STEPS (stdlib only)."""
import sys
src, dst, total = sys.argv[1], sys.argv[2], sys.argv[3]
L = open(src).read().split("\n"); L[0] = total
sp = [k for k, l in enumerate(L) if len(l.split()) == 20 and l.split()[0] == "1.000000"]; assert len(sp) == 1, sp
f = L[sp[0]].split(); f[10] = "0.070000"; L[sp[0]] = " ".join(f)
fl = [k for k, l in enumerate(L) if len(l.split()) == 32]; assert len(fl) == 1, fl
g = L[fl[0]].split(); g[4] = "100"; g[10] = "0"; g[11] = "100"; L[fl[0]] = " ".join(g)
open(dst, "w").write("\n".join(L)); print("switched", src, "->", dst, "steps", total)
