#!/usr/bin/env python3
"""Write emd_curves.dat for qhat_overlay.tex: the EMD qhat as the QUARK qhat.

The holographic qhat is the adjoint (gluon) one, qhat_A/T^3 = sqrt(lambda_t) Q(T)
(light-like adjoint Wilson loop, Liu-Rajagopal-Wiedemann); JETSCAPE quotes the
quark qhat, so qhat_q/T^3 = (C_F/C_A) sqrt(lambda_t) Q(T) = (4/9) sqrt(lambda_t) Q(T).
Q(T) is the fit of src/jet/EMDQhat.h to data/qhat_muB0.dat.
"""
import math


def Q(T):
    th = math.tanh((1000 * T - 133.6481) / 50.3029)
    return 0.5 * 3.98564 * (1 - th) + 0.5 * (1 + th) * (8.10975 + 4.10263 * math.exp(-(1000 * T - 133.6481) / 254.864))


cf = 4 / 9
rows = []
for i in range(71):
    T = 0.15 + 0.005 * i
    q = Q(T)
    rows.append((T, cf * q / 5, cf * q / 4, cf * q / 3, q / 4,
                 cf * math.sqrt(2.4) * q, cf * math.sqrt(2.0) * q, cf * math.sqrt(4.0) * q))
open("emd_curves.dat", "w").write(
    "# EMD q-hat is the adjoint (gluon) q-hat: q_A/T^3 = sqrt(lambda_t) Q(T) (Khan et al., LRW).\n"
    "# Quark q-hat = C_F/C_A q_A = (4/9) sqrt(lambda_t) Q(T), as plotted here against JETSCAPE (quark).\n"
    "# k25,k16,k9: lambda_t = 1/25, 1/16, 1/9 as quark; k16raw: lambda_t = 1/16 without C_F/C_A (as in Khan et al. Fig. 10)\n"
    "# ours: lambda_t = 2.4 (ourslo 2.0, ourshi 4.0) from the MATTER+LBT fit to ATLAS jets\n"
    "T k25 k16 k9 k16raw ours ourslo ourshi\n" + "".join(" ".join(f"{v:.4f}" for v in r) + "\n" for r in rows))
print("wrote emd_curves.dat")
