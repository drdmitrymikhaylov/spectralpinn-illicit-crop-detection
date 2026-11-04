"""Three things the 30 September check left open, on the same open spectra.

  1  operating points   the detector table quoted one point, 1 % false alarms.
                        Here the whole curve: detection at 0.1 / 0.5 / 1 / 2 /
                        5 / 10 % false alarms, and the false-alarm rate a survey
                        would have to accept for 90 / 95 / 99 % detection.
                        Plant-level bootstrap intervals on the 1 % column.
  2  honest threshold   the 1 % threshold in the table was set on the same
                        held-out leaves it was scored on.  Here it is set on
                        the training plants only (inner out-of-fold
                        posteriors), and the false-alarm rate actually
                        delivered on the test plants is read off.
  3  non-linear models  "the detector is a linear discriminant, a baseline".
                        A random forest and an RBF support-vector machine on
                        the same plant-grouped folds, same SNV spectra, same
                        400-1000 nm window: seven-way accuracy, the two grass
                        confusions, detection at 1 % false alarms.

Same data, same quality gate, same folds as species_separation_vnir.py (which
this script imports).  The LDA rows reproduce that script's numbers.

    python analysis/operating_points.py --data /path/to/data

Writes results/operating_points.json.  About a minute on a laptop.
"""

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.svm import SVC

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from species_separation_vnir import (FAR, SHORT, WINDOWS, auc, cv_posteriors,  # noqa: E402
                                     lda, physically_possible, snv)

FARS = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]
LEVELS = [0.90, 0.95, 0.99]
N_BOOT = 2000
SEED = 0


def load(data_dir):
    d = pathlib.Path(data_dir)
    meta = pd.read_csv(d / "crops_meta.csv", low_memory=False)
    z = np.load(d / "crops_spec.npz")
    wl, R = z["wl"].astype(float), z["R"].astype(float)
    keep = physically_possible(wl, R)
    meta, R = meta[keep].reset_index(drop=True), R[keep]
    species = meta["Common Name"].map(SHORT).to_numpy()
    names = sorted(set(species))
    y_all = np.array([names.index(s) for s in species])
    g = ((meta["Measurement Venue"] == "glasshouse")
         & (meta["Paired_Spectra"] == "leaf")).to_numpy()
    plant = (meta["Common Name"] + "|" + meta["Plant"].astype(str)).to_numpy()[g]
    return wl, R[g], y_all[g], plant, names


def detected_at(pos, neg, far):
    return float((pos > np.quantile(neg, 1 - far)).mean())


def operating_points(P, y, plant, names, rng):
    """Detection at several false-alarm rates, false alarms needed for several
    detection levels, and a plant-level bootstrap of the 1 % column, all from
    the out-of-fold posteriors."""
    plants = np.unique(plant)
    idx_by_plant = {p: np.flatnonzero(plant == p) for p in plants}
    out = {}
    for k, nm in enumerate(names):
        pos, neg = P[y == k, k], P[y != k, k]
        det = {str(f): detected_at(pos, neg, f) for f in FARS}
        need = {}
        for L in LEVELS:
            thr = np.quantile(pos, 1 - L)
            need[str(L)] = float((neg > thr).mean())
        boot = []
        for _ in range(N_BOOT):
            draw = rng.choice(plants, len(plants), replace=True)
            ii = np.concatenate([idx_by_plant[p] for p in draw])
            yy, pp = y[ii], P[ii, k]
            bpos, bneg = pp[yy == k], pp[yy != k]
            if len(bpos) == 0 or len(bneg) == 0:
                continue
            boot.append(detected_at(bpos, bneg, FAR))
        lo, hi = np.percentile(boot, [2.5, 97.5])
        out[nm] = {"n_target": int(len(pos)), "n_other": int(len(neg)),
                   "n_target_plants": int(len(np.unique(plant[y == k]))),
                   "auc": auc(pos, neg),
                   "detected_at_false_alarm_rate": det,
                   "false_alarm_rate_for_detection": need,
                   "detected_at_1pct_ci95_plants": [float(lo), float(hi)],
                   "n_bootstrap": len(boot)}
    return out


def honest_threshold(X, y, plant, names):
    """Threshold from the training plants (inner 5-fold out-of-fold posteriors,
    99th percentile of the non-target leaves), applied to the test plants."""
    counts = {nm: {"fa": 0, "n_other": 0, "det": 0, "n_target": 0} for nm in names}
    for tr, te in GroupKFold(n_splits=5).split(X, y, plant):
        P_in = cv_posteriors(X[tr], y[tr], plant[tr])
        P_te = lda(X[tr], y[tr]).predict_proba(X[te])
        for k, nm in enumerate(names):
            thr = np.quantile(P_in[y[tr] != k, k], 1 - FAR)
            c = counts[nm]
            c["fa"] += int((P_te[y[te] != k, k] > thr).sum())
            c["n_other"] += int((y[te] != k).sum())
            c["det"] += int((P_te[y[te] == k, k] > thr).sum())
            c["n_target"] += int((y[te] == k).sum())
    return {nm: {"false_alarm_rate_delivered": c["fa"] / c["n_other"],
                 "false_alarms": c["fa"], "n_other": c["n_other"],
                 "detected": c["det"] / c["n_target"],
                 "n_detected": c["det"], "n_target": c["n_target"]}
            for nm, c in counts.items()}


def oof_proba(model_fn, X, y, plant):
    P = np.zeros((len(y), len(np.unique(y))))
    for tr, te in GroupKFold(n_splits=5).split(X, y, plant):
        P[te] = model_fn().fit(X[tr], y[tr]).predict_proba(X[te])
    return P


def seven_way(P, y, names):
    pred = P.argmax(1)
    cm = np.zeros((len(names), len(names)), int)
    for t, p in zip(y, pred):
        cm[t, p] += 1
    i = {nm: k for k, nm in enumerate(names)}
    det = {nm: {"auc": auc(P[y == k, k], P[y != k, k]),
                "detected_at_1pct_false_alarm":
                    detected_at(P[y == k, k], P[y != k, k], FAR)}
           for k, nm in enumerate(names)}
    return {"accuracy": float((pred == y).mean()), "n_errors": int((pred != y).sum()),
            "millet_called_sorghum": int(cm[i["millet"], i["sorghum"]]),
            "sorghum_called_millet": int(cm[i["sorghum"], i["millet"]]),
            "poplar_called_pumpkin": int(cm[i["poplar"], i["pumpkin"]]),
            "recall": {nm: float(cm[k, k] / cm[k].sum()) for k, nm in enumerate(names)},
            "detector": det, "confusion_rows_true": cm.tolist()}


MODELS = {
    "random forest": lambda: RandomForestClassifier(n_estimators=500, random_state=SEED,
                                                    n_jobs=-1),
    "rbf svm": lambda: SVC(kernel="rbf", C=10.0, gamma="scale", probability=True,
                           random_state=SEED),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    a = ap.parse_args()
    wl, R, y, plant, names = load(a.data)
    rng = np.random.default_rng(SEED)

    out = {"n_leaves": int(len(y)), "n_plants": int(len(np.unique(plant))),
           "false_alarm_rates": FARS, "detection_levels": LEVELS,
           "n_bootstrap": N_BOOT, "operating_points": {}, "honest_threshold": {},
           "models": {}}
    for wname, (lo, hi) in WINDOWS.items():
        w = (wl >= lo) & (wl <= hi)
        X = snv(R[:, w])
        P = cv_posteriors(X, y, plant)
        out["operating_points"][wname] = operating_points(P, y, plant, names, rng)
        out["honest_threshold"][wname] = honest_threshold(X, y, plant, names)

    w = (wl >= 400) & (wl <= 1000)
    X = snv(R[:, w])
    out["models"]["window"] = "400-1000 nm"
    out["models"]["linear discriminant"] = seven_way(cv_posteriors(X, y, plant), y, names)
    for mname in ("random forest", "rbf svm"):
        out["models"][mname] = seven_way(oof_proba(MODELS[mname], X, y, plant), y, names)
    out["models"]["settings"] = {
        "random forest": "500 trees, scikit-learn defaults otherwise, seed 0",
        "rbf svm": "C = 10, gamma = 'scale', Platt probabilities, seed 0",
        "note": "untuned; the same plant-grouped 5 folds and SNV spectra as the LDA"}

    pathlib.Path("results").mkdir(exist_ok=True)
    with open("results/operating_points.json", "w") as fh:
        json.dump(out, fh, indent=1)
    brief = json.loads(json.dumps(out))
    for m in ("linear discriminant", "random forest", "rbf svm"):
        brief["models"][m].pop("confusion_rows_true")
    print(json.dumps(brief, indent=1))


if __name__ == "__main__":
    main()
