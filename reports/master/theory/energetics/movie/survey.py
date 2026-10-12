import sys, numpy as np
sys.argv = [sys.argv[0], sys.argv[1], "0"]
exec(open("escape_movie.py").read().split('print(len(ev), "depleting escapes")')[0])
for q, e in enumerate(ev):
    ch = follow(e); n = e["n"]; full = all(m in ch for m in range(n - 30, n + 13))
    if not full: print(q, "incomplete"); continue
    c = [cost(ch[m], m)[0] for m in (n - 30, n - 20, n - 10, n)]; P = np.array([ch[m] for m in range(n - 30, n + 1)])
    jit = np.hypot(*np.diff(P, axis=0).T).sum(); dep = tree.query(ch[n + 12])[0]
    print(q, "cos", np.round(c, 2), "path", round(jit, 1), "net", round(np.hypot(*(P[-1] - P[0])), 1), "depth+24", round(dep, 1))
