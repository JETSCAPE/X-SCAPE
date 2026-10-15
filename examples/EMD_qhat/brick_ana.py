#!/usr/bin/env python3
"""Parton-level brick observables: leading-parton pT and cone (R) pT around the
gun axis (eta=0), counting status>=0 partons and subtracting status=-1 holes."""
import math, sys

def events(path):
    ev = []
    for line in open(path):
        if line.startswith('#'):
            if 'Event' in line and ev is not None:
                if ev: yield ev
                ev = []
            continue
        f = line.split()
        if len(f) < 7: continue
        pid, st = int(f[1]), int(f[2]); E, px, py, pz = map(float, f[3:7])
        ev.append((pid, st, E, px, py, pz))
    if ev: yield ev

def kin(px, py, pz):
    pt = math.hypot(px, py); p = math.sqrt(pt*pt + pz*pz)
    eta = 0.5*math.log((p+pz)/(p-pz)) if p > abs(pz) else 0.
    return pt, eta, math.atan2(py, px)

def analyse(path, R=0.4):
    lead, cone, ngl, n = 0., 0., 0, 0
    for ev in events(path):
        n += 1
        parts = [(st, *kin(px, py, pz)) for pid, st, E, px, py, pz in ev]
        pos = [p for p in parts if p[0] >= 0]
        l = max(pos, key=lambda p: p[1]); lead += l[1]
        c = 0.
        for st, pt, eta, phi in parts:
            dphi = abs((phi - l[3] + math.pi) % (2*math.pi) - math.pi)
            if math.hypot(eta - l[2], dphi) < R:
                c += pt if st >= 0 else (-pt if st == -1 else 0.)  # -13: momentum given to the medium
        cone += c; ngl += len(pos)
    return n, lead/n, cone/n, ngl/n

if __name__ == '__main__':
    print(f"{'run':14s} {'nev':>5s} {'<lead pT>':>10s} {'<cone pT R=0.4>':>16s} {'<N partons>':>12s}")
    for f in sys.argv[1:]:
        n, l, c, g = analyse(f)
        print(f"{f.split('/')[-1].replace('_final_state_partons.dat',''):14s} {n:5d} {l:10.2f} {c:16.2f} {g:12.1f}")
