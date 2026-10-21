#!/usr/bin/env python3
"""Check of frank_constant.py (2026-10-21): ML convergence from two starting points; spectrum plateau at |q| > 2.5 (bound on per-cell director noise);
ML fit with an additive white-noise term, mean P = A/(K1 qperp^2 + K3 qpar^2) + A sigma^2, at qmax 0.5, 0.8, 1.2. Usage: frank_check.py RUN"""
import sys, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/Kcalib")
import frank_constant as F
from scipy.optimize import minimize
run = sys.argv[1]
frames = F.read_frames(run, 200.0); L, th0, P, S = F.spectra(frames); nf = len(frames); A = L * L
q = 2*np.pi*np.fft.fftfreq(L); QX, QY = np.meshgrid(q, q, indexing="ij"); Q2 = QX**2+QY**2
out = []
# per-cell director noise estimate from the spectrum plateau at large q (|q| > 2.5)
hi = Q2 > 2.5**2
out.append(f"{run.split('/')[-1]}: plateau <P>/A at |q|>2.5 = {P[:, hi].mean()/A:.4f} (white-noise sigma^2 upper bound); expected 1/(K q^2) at q=0.8: bend {1/(22*0.64):.3f}, splay {1/(43*0.64):.3f}")
for qmax in (0.5, 0.8, 1.2):
    m = (Q2 > 0) & (Q2 < qmax**2); qx, qy = QX[m], QY[m]; Pm = P[:, m]
    qp, qn = F.rotate(qx[None], qy[None], th0[:, None]); a, b, Pf = (qp**2).ravel(), (qn**2).ravel(), Pm.ravel()
    # convergence: two starting points, report success
    def nll(p):
        lam = (p[0]*b + p[1]*a + p[2]) / A; r = -1/lam + Pf
        return np.sum(-np.log(lam) + lam*Pf), np.array([np.sum(r*b), np.sum(r*a), np.sum(r)])/A
    res = [minimize(nll, x0=x0, jac=True, method="L-BFGS-B", bounds=[(1e-3,None)]*2+[(0,None)], options=dict(ftol=1e-15, gtol=1e-10, maxiter=5000)) for x0 in ([30,30,0.05],[60,10,0.5])]
    # noise model: mean P = A/(K1 qn^2 + K3 qp^2) + A*s2
    def nll2(p):
        mu = A/(p[0]*b + p[1]*a) + A*p[2]; return np.sum(np.log(mu) + Pf/mu)
    r2 = minimize(nll2, x0=[43, 22, 0.003], method="L-BFGS-B", bounds=[(1,None),(1,None),(0,None)], options=dict(ftol=1e-15, maxiter=5000))
    out.append(f"  qmax {qmax}: dir_ml x0a K1 {res[0].x[0]:.2f} K3 {res[0].x[1]:.2f} c {res[0].x[2]:.3g} ok={res[0].success} | x0b K1 {res[1].x[0]:.2f} K3 {res[1].x[1]:.2f} ok={res[1].success} | noise model K1 {r2.x[0]:.2f} K3 {r2.x[1]:.2f} ratio {r2.x[0]/r2.x[1]:.2f} sigma2 {r2.x[2]:.4f}")
print("\n".join(out))
