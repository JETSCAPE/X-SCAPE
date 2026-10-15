"""Read and combine jetana spectrum files; R_AA with statistical errors."""
import numpy as np


def read_spec(path):
    d = {}
    for line in open(path):
        f = line.split()
        if f[0] == "sumw":
            d["sumw"] = float(f[1])
        else:
            d[f[0]] = np.array([float(x) for x in f[1:]])
    return d


def combine(paths):
    tot = None
    for p in paths:
        d = read_spec(p)
        if tot is None:
            tot = d
            continue
        for k in d:
            if k.endswith("_h") or k.endswith("_h2") or k == "sumw":
                tot[k] = tot[k] + d[k]
    return tot


def raa(aa, pp, obs):
    ya, yp = aa[obs + "_h"] / aa["sumw"], pp[obs + "_h"] / pp["sumw"]
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ya / yp
        e = r * np.sqrt(aa[obs + "_h2"] / aa[obs + "_h"] ** 2 + pp[obs + "_h2"] / pp[obs + "_h"] ** 2)
    return aa[obs], r, e
