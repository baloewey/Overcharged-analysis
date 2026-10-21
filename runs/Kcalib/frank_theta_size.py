#!/usr/bin/env python3
"""Small-angle check for frank_constant.py (2026-10-21): per-cell and long-wavelength (|q| < 0.8) rms of the director deviation, kurtosis, and the
power of the terms dropped in the linearisation (-sin(theta) d_x theta in div n, sin(theta) d_y theta in curl n) relative to the linear ones. Usage: frank_theta_size.py RUN"""
import sys, numpy as np
sys.path.insert(0, "/home/baloewe/Physics/Overcharged/runs/Kcalib")
import frank_constant as F
run = sys.argv[1]; fr = F.read_frames(run, 200.0)
L = int(fr[0][:,0].max())+1; A=L*L
d2, d4, big, sm2, nl = [], [], [], [], []
q = 2*np.pi*np.fft.fftfreq(L); QX,QY=np.meshgrid(q,q,indexing="ij"); keep = (QX**2+QY**2) < 0.8**2
for f in fr[::5]:
    ix, iy = f[:,0].astype(int), f[:,1].astype(int); ang = np.arctan2(f[:,3], f[:,2])
    th0 = 0.5*np.arctan2(np.mean(np.sin(2*ang)), np.mean(np.cos(2*ang)))
    d = (ang-th0+np.pi/2)%np.pi-np.pi/2; th=np.zeros((L,L)); th[ix,iy]=d
    d2.append(np.mean(d**2)); d4.append(np.mean(d**4)); big.append(np.mean(np.abs(d)>np.radians(30)))
    # long-wavelength part only (|q| < 0.8, the modes the fit uses)
    ths = np.real(np.fft.ifft2(np.fft.fft2(th)*keep)); sm2.append(np.mean(ths**2))
    # exact vs linear splay/bend from the smoothed field, director frame (x along th0 after rotation is approximated: th0 small, use lab grid)
    c, s = np.cos(ths), np.sin(ths)
    dx = lambda a: (np.roll(a,-1,0)-np.roll(a,1,0))/2; dy = lambda a: (np.roll(a,-1,1)-np.roll(a,1,1))/2
    div_lin, curl_lin = dy(ths), dx(ths); div_ex = -s*dx(ths)+c*dy(ths); curl_ex = c*dx(ths)+s*dy(ths)
    nl.append((np.mean((div_ex-div_lin)**2)/np.mean(div_lin**2), np.mean((curl_ex-curl_lin)**2)/np.mean(curl_lin**2),
               np.mean(div_ex**2)/np.mean(div_lin**2)-1, np.mean(curl_ex**2)/np.mean(curl_lin**2)-1))
nl = np.array(nl).mean(0)
print(f"{run.split('/')[-1]}: per-cell theta rms {np.degrees(np.sqrt(np.mean(d2))):.1f} deg, <theta^4>/<theta^2>^2 {np.mean(d4)/np.mean(d2)**2:.2f}, frac |theta|>30deg {np.mean(big):.4f}; "
      f"long-wavelength (|q|<0.8) theta rms {np.degrees(np.sqrt(np.mean(sm2))):.2f} deg, <sin^2> {np.mean(sm2):.4f}; "
      f"neglected/linear power: splay {nl[0]:.4f}, bend {nl[1]:.4f}; exact/linear -1: splay {nl[2]:+.4f}, bend {nl[3]:+.4f}")
