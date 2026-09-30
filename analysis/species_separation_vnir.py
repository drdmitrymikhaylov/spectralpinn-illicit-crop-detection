"""How well does the 400-1000 nm window separate one crop species from another,
and how much of a pixel must the target fill?

A reference calculation on OPEN leaf spectra.  It uses no survey data and none
of the species is a target of the field work this repository describes; the
seven crops stand in for "a species among other species".

Data: A. C. Burnett, S. P. Serbin, K. J. Davidson, K. S. Ely, A. Rogers (2020),
"Hyperspectral leaf reflectance, biochemistry, and physiology of droughted and
watered crops", EcoSIS package 5090905b-176c-4d17-bf60-59a69939eea6.  One
spectroradiometer, seven species, leaf clip.  Not redistributed here: download
from https://ecosis.org, resample reflectance to 4 nm over 400-2400 nm and
save as crops_spec.npz (arrays `wl`, `R`) next to the metadata table
crops_meta.csv.

Four questions.

  1  separation      seven-way species call from glasshouse leaves, test plants
                     never seen in training; again with whole measurement
                     dates held out; 400-1000 nm against 400-2400 nm
  2  detector        each species in turn as the target: detection rate at a
                     1 % false-alarm rate against the other six
  3  mixed pixel     target and host leaf mixed linearly by area fraction f:
                     detection at 1 % false alarms against f, all 42 ordered
                     species pairs
  4  new venue       the glasshouse model asked about the same species (pumpkin)
                     measured in a field, leaf clip and above-canopy

    python analysis/species_separation_vnir.py --data /path/to/data

Writes results/species_separation_vnir.json.  About 20 s on a laptop.
"""

import argparse
import json
import pathlib

import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import GroupKFold

WINDOWS = {"400-1000 nm": (400, 1000), "400-2400 nm": (400, 2400)}
FRACTIONS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
FAR = 0.01
N_MIX = 400
SHORT = {"Carolina poplar": "poplar", "cayenne pepper": "pepper",
         "common sunflower": "sunflower", "cultivated radish": "radish",
         "field pumpkin": "pumpkin", "foxtail millet": "millet", "sorghum": "sorghum"}


def physically_possible(wl, R):
    """Reflectance inside [0, 1], a near-infrared plateau above the red trough,
    no detector-splice step between adjacent bands."""
    ok = np.isfinite(R).all(1) & (R.min(1) >= 0) & (R.max(1) <= 1)
    red = R[:, (wl >= 660) & (wl <= 680)].mean(1)
    nir = R[:, (wl >= 800) & (wl <= 900)].mean(1)
    return ok & (nir > red) & (np.abs(np.diff(R, axis=1)).max(1) <= 0.25)


def snv(X):
    """Per-spectrum standardisation: removes brightness, keeps shape."""
    return (X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True)


def lda(X, y):
    return LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto").fit(X, y)


def auc(pos, neg):
    allv = np.concatenate([pos, neg])
    r = np.argsort(np.argsort(allv)) + 1.0
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2)
                 / (len(pos) * len(neg)))


def cv_posteriors(X, y, groups, n_splits=5):
    P = np.zeros((len(y), len(np.unique(y))))
    for tr, te in GroupKFold(n_splits=n_splits).split(X, y, groups):
        P[te] = lda(X[tr], y[tr]).predict_proba(X[te])
    return P


def separation(X, y, groups, names):
    P = cv_posteriors(X, y, groups)
    pred = P.argmax(1)
    cm = np.zeros((len(names), len(names)), int)
    for t, p in zip(y, pred):
        cm[t, p] += 1
    recall = cm.diagonal() / cm.sum(1)
    out = {"n": int(len(y)), "accuracy": float((pred == y).mean()),
           "balanced_accuracy": float(recall.mean()),
           "recall": {names[k]: float(recall[k]) for k in range(len(names))},
           "confusion_rows_true": cm.tolist(), "n_errors": int((pred != y).sum())}
    return out, P


def detector(P, y, names):
    """One species against the other six, from the out-of-fold posteriors."""
    out = {}
    for k, nm in enumerate(names):
        pos, neg = P[y == k, k], P[y != k, k]
        thr = np.quantile(neg, 1 - FAR)
        out[nm] = {"n_target": int(len(pos)), "n_other": int(len(neg)),
                   "auc": auc(pos, neg),
                   "detected_at_1pct_false_alarm": float((pos > thr).mean()),
                   "missed_of_target": int((pos <= thr).sum())}
    return out


def mixed_pixel(Rw, y, groups, names, seed=0):
    """Target leaf and host leaf mixed by area fraction f in reflectance, then
    standardised as one spectrum.  The detector is trained on pure leaves of
    the two species from other plants; the threshold gives 1 % false alarms on
    pure host leaves of the test plants."""
    rng = np.random.default_rng(seed)
    pairs = {}
    for t in range(len(names)):
        for h in range(len(names)):
            if t == h:
                continue
            m = (y == t) | (y == h)
            R2, y2, g2 = Rw[m], (y[m] == t).astype(int), groups[m]
            hit = {f: [] for f in FRACTIONS}
            for tr, te in GroupKFold(n_splits=5).split(R2, y2, g2):
                model = lda(snv(R2[tr]), y2[tr])
                T, H = R2[te][y2[te] == 1], R2[te][y2[te] == 0]
                if len(T) == 0 or len(H) == 0:
                    continue
                thr = np.quantile(model.decision_function(snv(H)), 1 - FAR)
                i, j = rng.integers(0, len(T), N_MIX), rng.integers(0, len(H), N_MIX)
                for f in FRACTIONS:
                    mix = f * T[i] + (1 - f) * H[j]
                    hit[f].append(float((model.decision_function(snv(mix)) > thr).mean()))
            rate = {f: float(np.mean(v)) for f, v in hit.items()}
            pairs[f"{names[t]} in {names[h]}"] = {
                "detected": {str(f): rate[f] for f in FRACTIONS},
                "fraction_for_90pct": fraction_for(rate, 0.9),
                "fraction_for_50pct": fraction_for(rate, 0.5)}
    summary = {}
    for f in FRACTIONS:
        v = np.array([p["detected"][str(f)] for p in pairs.values()])
        summary[str(f)] = {"median": float(np.median(v)), "min": float(v.min()),
                           "max": float(v.max()),
                           "pairs_at_or_above_90pct": int((v >= 0.9).sum())}
    for key in ("fraction_for_90pct", "fraction_for_50pct"):
        v = np.array([p[key] if p[key] is not None else np.nan for p in pairs.values()])
        summary[key] = {"median": float(np.nanmedian(v)), "min": float(np.nanmin(v)),
                        "max": float(np.nanmax(v)),
                        "pairs_never_reaching_it": int(np.isnan(v).sum())}
    worst = min(pairs, key=lambda k: pairs[k]["detected"]["0.5"])
    best = max(pairs, key=lambda k: pairs[k]["detected"]["0.2"])
    return {"n_pairs": len(pairs), "n_mixtures_per_fold": N_MIX,
            "false_alarm_rate": FAR, "fractions": FRACTIONS, "summary": summary,
            "hardest_pair_at_half": worst, "easiest_pair_at_fifth": best,
            "pairs": pairs}


def fraction_for(rate, level):
    """Smallest area fraction reaching `level`, by linear interpolation."""
    fs = sorted(rate)
    prev_f, prev_r = 0.0, FAR          # at f = 0 the pixel is pure host
    for f in fs:
        if rate[f] >= level:
            if rate[f] == prev_r:
                return float(f)
            return float(prev_f + (level - prev_r) * (f - prev_f) / (rate[f] - prev_r))
        prev_f, prev_r = f, rate[f]
    return None


def new_venue(X_train, y_train, X_field, names, target):
    model = lda(X_train, y_train)
    pred = model.predict(X_field)
    k = names.index(target)
    called = {names[c]: int((pred == c).sum()) for c in range(len(names))}
    return {"n": int(len(pred)), "called_target": int((pred == k).sum()),
            "rate": float((pred == k).mean()), "called": called}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    a = ap.parse_args()
    d = pathlib.Path(a.data)
    meta = pd.read_csv(d / "crops_meta.csv", low_memory=False)
    z = np.load(d / "crops_spec.npz")
    wl, R = z["wl"].astype(float), z["R"].astype(float)
    keep = physically_possible(wl, R)
    meta, R = meta[keep].reset_index(drop=True), R[keep]

    species = meta["Common Name"].map(SHORT).to_numpy()
    names = sorted(set(species))
    y_all = np.array([names.index(s) for s in species])
    glass = (meta["Measurement Venue"] == "glasshouse").to_numpy()
    leaf = (meta["Paired_Spectra"] == "leaf").to_numpy()
    plant = (meta["Common Name"] + "|" + meta["Plant"].astype(str)).to_numpy()
    date = meta["Measurement Date"].astype(str).to_numpy()

    g = glass & leaf
    y = y_all[g]
    out = {"data": {"n_after_quality_gate": int(len(R)),
                    "n_dropped_by_quality_gate": int((~keep).sum()),
                    "n_glasshouse_leaf": int(g.sum()),
                    "per_species": {nm: int((y == k).sum()) for k, nm in enumerate(names)},
                    "plants": int(len(np.unique(plant[g]))),
                    "dates": int(len(np.unique(date[g]))),
                    "band_step_nm": float(np.median(np.diff(wl)))},
           "separation": {}, "detector": {}}

    for wname, (lo, hi) in WINDOWS.items():
        w = (wl >= lo) & (wl <= hi)
        X = snv(R[g][:, w])
        by_plant, P = separation(X, y, plant[g], names)
        by_date, _ = separation(X, y, date[g], names)
        out["separation"][wname] = {"n_bands": int(w.sum()),
                                    "test_plants_unseen": by_plant,
                                    "test_dates_unseen": by_date}
        out["detector"][wname] = detector(P, y, names)

    w = (wl >= 400) & (wl <= 1000)
    out["mixed_pixel"] = mixed_pixel(R[g][:, w], y, plant[g], names)

    X = snv(R[g][:, w])
    field = ~glass
    out["new_venue"] = {
        "target": "pumpkin", "window": "400-1000 nm",
        "glasshouse_recall_test_plants_unseen":
            out["separation"]["400-1000 nm"]["test_plants_unseen"]["recall"]["pumpkin"],
        "field_leaf": new_venue(X, y, snv(R[field & leaf][:, w]), names, "pumpkin"),
        "field_canopy": new_venue(X, y, snv(R[field & ~leaf][:, w]), names, "pumpkin")}

    pathlib.Path("results").mkdir(exist_ok=True)
    with open("results/species_separation_vnir.json", "w") as fh:
        json.dump(out, fh, indent=1)
    brief = json.loads(json.dumps(out))
    brief["mixed_pixel"].pop("pairs")
    print(json.dumps(brief, indent=1))


if __name__ == "__main__":
    main()
