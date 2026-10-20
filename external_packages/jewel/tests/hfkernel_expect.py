"""hfkernel_expect.py <test_hfkernel output>: exact expectation for test_hfkernel by numerical integration.

Emission density of an eikonal source = (dt/t) dz alpha_s/(2pi) P(z,t) over the same phase space as getmass:
  light: Q0^2 <= t <= E^2, z in [eps, 1-eps], eps = 1/2 - 1/2 sqrt(1-Q0^2/t) sqrt(1-t/E^2), P = C_F (1+z^2)/(1-z),
         alpha_s(z(1-z) t)
  heavy: Q0^2 <= t <= E^2 - m^2, z in two-body limits with daughter masses sqrt(m^2+Q0^2/4), Q0/2,
         P = C_F [(1+z^2)/(1-z) - 2 m^2/t], alpha_s(max(kt^2, 3 Q0^2/16)), kt^2 = z(1-z)t - (1-z)^2 m^2
histogrammed in ln(1/theta), theta = kt/(z(1-z)E); compared with the MC columns.
"""
import sys

import numpy as np

f = sys.argv[1]
hdr = open(f).readline().split()
E = float(hdr[1].split("=")[1])
m = float(hdr[2].split("=")[1])
mc = np.loadtxt(f)
Q0, LPS, CF = 1.5, 0.4, 4.0 / 3.0


def alphas(t):
    return 4 * np.pi / ((11 - 2 * 3 / 3) * np.log(t / LPS**2))


def zlim(t, mq, m1, m2):
    M2 = t + mq * mq
    lam = (1 - (m1 + m2) ** 2 / M2) * (1 - (m1 - m2) ** 2 / M2)
    beta = np.sqrt(np.clip(1 - M2 / E**2, 0, None))
    c = 0.5 * (1 + (m1 * m1 - m2 * m2) / M2)
    h = 0.5 * beta * np.sqrt(np.clip(lam, 0, None))
    ok = (lam > 0) & (M2 < E**2)
    return c - h, c + h, ok


def density(mq, nt=4000, nz=4000):
    edges = np.linspace(0, 5, 41)
    hist = np.zeros(40)
    tmax = E * E - mq * mq
    lt = np.linspace(np.log(Q0**2), np.log(tmax), nt + 1)
    tc = np.exp(0.5 * (lt[1:] + lt[:-1]))
    dlt = lt[1] - lt[0]
    for t in tc:
        if mq > 0:
            zlo, zhi, ok = zlim(t, mq, np.sqrt(mq * mq + Q0**2 / 4), Q0 / 2)
        else:
            eps = 0.5 - 0.5 * np.sqrt(1 - Q0**2 / t) * np.sqrt(1 - t / E**2)
            zlo, zhi, ok = eps, 1 - eps, t < E * E
        if not ok or zhi <= zlo:
            continue
        # integrate in u = -ln(1-z) (the kernel is ~ 1/(1-z))
        ulo, uhi = -np.log(1 - zlo), -np.log(1 - zhi)
        u = np.linspace(ulo, uhi, nz + 1)
        uc = 0.5 * (u[1:] + u[:-1])
        z = 1 - np.exp(-uc)
        du = u[1] - u[0]
        jac = 1 - z  # dz = (1-z) du
        if mq > 0:
            kt2 = z * (1 - z) * t - (1 - z) ** 2 * mq * mq
            P = CF * ((1 + z * z) / (1 - z) - 2 * mq * mq / t)
            a = alphas(np.maximum(kt2, 3 * Q0**2 / 16))
        else:
            kt2 = z * (1 - z) * t
            P = CF * (1 + z * z) / (1 - z)
            a = alphas(z * (1 - z) * t)
        w = a / (2 * np.pi) * np.clip(P, 0, None) * jac * du * dlt  # dt/t = dln t
        good = kt2 > 0
        th = np.sqrt(kt2[good]) / (z[good] * (1 - z[good]) * E)
        hist += np.histogram(np.log(1 / th), bins=edges, weights=w[good])[0]
    return hist / (edges[1] - edges[0])


dl, dh = density(0.0), density(m)
print(f"# E={E} m={m}: ln(1/theta)  light MC  light exact | heavy MC  heavy exact | ratio MC  ratio exact")
for i, row in enumerate(mc):
    if row[1] == 0 and dl[i] == 0:
        continue
    rx = dh[i] / dl[i] if dl[i] > 0 else float("nan")
    print(f"{row[0]:6.3f}   {row[1]:7.4f} {dl[i]:7.4f}   | {row[2]:7.4f} {dh[i]:7.4f}   | {row[3]:6.3f}+-{row[4]:5.3f} {rx:6.3f}")
