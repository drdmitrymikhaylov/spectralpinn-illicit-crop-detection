"""Pins the README section "Operating points, a threshold set without the test
leaves, and two non-linear models" to results/operating_points.json, and that
file's linear-discriminant rows to results/species_separation_vnir.json.
No dataset needed."""

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
R = json.loads((ROOT / "results" / "operating_points.json").read_text())
S = json.loads((ROOT / "results" / "species_separation_vnir.json").read_text())
README = " ".join((ROOT / "README.md").read_text(encoding="utf-8").split())
NAMES = ["pepper", "radish", "sunflower", "poplar", "pumpkin", "sorghum", "millet"]
W = "400-1000 nm"
FARS = ["0.001", "0.005", "0.01", "0.02", "0.05", "0.1"]


def pct(v):
    return f"{100 * v:.1f}%"


def cell(v):
    return v if abs(v) >= 1e-12 else 0.0


def test_lda_rows_reproduce_the_september_run():
    for w in (W, "400-2400 nm"):
        for nm in NAMES:
            a, b = R["operating_points"][w][nm], S["detector"][w][nm]
            assert abs(a["auc"] - b["auc"]) < 1e-9
            assert abs(a["detected_at_false_alarm_rate"]["0.01"]
                       - b["detected_at_1pct_false_alarm"]) < 1e-9
            assert a["n_target"] == b["n_target"] and a["n_other"] == b["n_other"]
    m = R["models"]["linear discriminant"]
    s = S["separation"][W]["test_plants_unseen"]
    assert abs(m["accuracy"] - s["accuracy"]) < 1e-9
    assert m["n_errors"] == s["n_errors"] == 168
    assert m["confusion_rows_true"] == s["confusion_rows_true"]


def test_operating_curve_table():
    op = R["operating_points"][W]
    for nm in NAMES:
        o = op[nm]
        d, ci = o["detected_at_false_alarm_rate"], o["detected_at_1pct_ci95_plants"]
        row = (f"| {nm} | " + " | ".join(pct(d[k]) for k in FARS)
               + f" | [{pct(ci[0])}, {pct(ci[1])}] |")
        assert row in README, row
        # detection is monotone in the false-alarm rate
        vals = [d[k] for k in FARS]
        assert all(x <= y + 1e-12 for x, y in zip(vals, vals[1:]))
        assert ci[0] <= d["0.01"] <= ci[1]
        assert o["n_bootstrap"] == R["n_bootstrap"] == 2000
    assert R["n_plants"] == 128 and R["n_leaves"] == 2009
    assert op["millet"]["n_target_plants"] == 16 and op["pumpkin"]["n_target_plants"] == 24


def test_operating_curve_prose():
    op = R["operating_points"][W]
    d = {nm: op[nm]["detected_at_false_alarm_rate"] for nm in NAMES}
    assert pct(d["millet"]["0.01"]) == "53.1%" and pct(d["millet"]["0.02"]) == "82.3%"
    assert pct(d["millet"]["0.05"]) == "97.9%"
    assert min(d[nm]["0.05"] for nm in NAMES) > 0.94
    assert min(d[nm]["0.1"] for nm in NAMES) > 0.98
    assert pct(d["pepper"]["0.001"]) == "98.8%" and pct(d["millet"]["0.001"]) == "13.5%"
    assert pct(d["sunflower"]["0.005"]) == "94.2%" and pct(d["sunflower"]["0.001"]) == "62.2%"
    ci = {nm: op[nm]["detected_at_1pct_ci95_plants"] for nm in NAMES}
    assert pct(ci["millet"][0]) == "32.4%" and pct(ci["millet"][1]) == "90.7%"
    assert 0.30 < ci["millet"][0] < 0.34 and 0.89 < ci["millet"][1] < 0.91
    assert pct(ci["sorghum"][0]) == "44.1%" and pct(ci["sorghum"][1]) == "90.4%"
    assert all(ci[nm][0] > 0.85 for nm in ("pepper", "radish", "sunflower", "poplar"))
    assert pct(ci["pumpkin"][0]) == "80.1%" and pct(ci["pumpkin"][1]) == "93.3%"
    assert "resampling plants puts it anywhere from a third to nine in ten" in README


def test_cost_and_honest_threshold_table():
    op, ht = R["operating_points"][W], R["honest_threshold"][W]
    for nm in NAMES:
        n, h = op[nm]["false_alarm_rate_for_detection"], ht[nm]
        row = (f"| {nm} | {pct(n['0.9'])} | {pct(n['0.95'])} | {pct(n['0.99'])} | "
               f"{pct(h['false_alarm_rate_delivered'])} ({h['false_alarms']} of "
               f"{h['n_other']:,}".replace(",", " ") + f") | {pct(h['detected'])} |")
        assert row in README, row
        assert h["false_alarms"] / h["n_other"] == pytest.approx(h["false_alarm_rate_delivered"])
        assert h["n_detected"] / h["n_target"] == pytest.approx(h["detected"])
        assert h["n_target"] == op[nm]["n_target"] and h["n_other"] == op[nm]["n_other"]
        assert n["0.9"] <= n["0.95"] <= n["0.99"]


def test_cost_and_honest_threshold_prose():
    op, ht = R["operating_points"][W], R["honest_threshold"][W]
    n = {nm: op[nm]["false_alarm_rate_for_detection"] for nm in NAMES}
    assert pct(n["millet"]["0.9"]) == "3.0%" and pct(n["sorghum"]["0.9"]) == "2.9%"
    assert pct(n["pumpkin"]["0.9"]) == "1.8%"
    assert all(round(100 * n[nm]["0.9"], 1) <= 0.6
               for nm in ("pepper", "radish", "sunflower", "poplar"))
    assert pct(n["pumpkin"]["0.99"]) == "18.7%"
    assert pct(n["millet"]["0.99"]) == "5.2%" and pct(n["sorghum"]["0.99"]) == "5.3%"
    assert all(round(100 * n[nm]["0.99"], 1) <= 2.6
               for nm in ("pepper", "radish", "sunflower", "poplar"))
    fa = [ht[nm]["false_alarm_rate_delivered"] for nm in NAMES]
    assert pct(min(fa)) == "0.5%" and pct(max(fa)) == "1.3%"
    assert max(fa) == ht["pumpkin"]["false_alarm_rate_delivered"]
    assert ht["pumpkin"]["false_alarms"] == 20 and ht["pumpkin"]["n_other"] == 1542
    moves = {nm: ht[nm]["detected"] - op[nm]["detected_at_false_alarm_rate"]["0.01"]
             for nm in NAMES}
    assert max(abs(v) for v in moves.values()) == pytest.approx(moves["millet"])
    assert round(100 * moves["millet"], 1) == 4.2
    assert pct(ht["millet"]["detected"]) == "57.3%" and pct(ht["poplar"]["detected"]) == "92.0%"
    assert pct(ht["sorghum"]["detected"]) == "73.5%" and pct(ht["radish"]["detected"]) == "96.3%"
    for nm in ("pepper", "sunflower", "pumpkin"):
        assert pct(ht[nm]["detected"]) == pct(op[nm]["detected_at_false_alarm_rate"]["0.01"])
    full = R["honest_threshold"]["400-2400 nm"]
    fa2 = [full[nm]["false_alarm_rate_delivered"] for nm in NAMES]
    det2 = [full[nm]["detected"] for nm in NAMES]
    assert pct(min(fa2)) == "0.6%" and pct(max(fa2)) == "1.1%"
    assert pct(min(det2)) == "98.9%" and pct(max(det2)) == "100.0%"
    assert full["millet"]["n_detected"] == 95 and full["millet"]["n_target"] == 96


def test_models_table_and_prose():
    M = R["models"]
    assert M["window"] == W
    rows = {"linear discriminant": "linear discriminant", "random forest": "random forest",
            "rbf svm": "RBF SVM"}
    for key, label in rows.items():
        m = M[key]
        grass = m["millet_called_sorghum"] + m["sorghum_called_millet"]
        det = m["detector"]
        row = (f"| {label} | {pct(m['accuracy'])} | {m['n_errors']} | {grass} | "
               f"{m['poplar_called_pumpkin']} | "
               + " | ".join(pct(det[nm]["detected_at_1pct_false_alarm"])
                            for nm in ("millet", "sorghum", "pumpkin", "poplar")) + " |")
        assert row in README, row
        assert m["n_errors"] == 2009 - sum(m["confusion_rows_true"][k][k] for k in range(7))
    lda, rf, svm = M["linear discriminant"], M["random forest"], M["rbf svm"]
    assert lda["n_errors"] - rf["n_errors"] == 32
    assert lda["poplar_called_pumpkin"] - rf["poplar_called_pumpkin"] == 26
    assert rf["millet_called_sorghum"] + rf["sorghum_called_millet"] == 42
    assert lda["millet_called_sorghum"] + lda["sorghum_called_millet"] == 55
    assert pct(rf["detector"]["millet"]["detected_at_1pct_false_alarm"]) == "78.1%"
    for nm in ("poplar", "pumpkin", "radish"):
        assert (rf["detector"][nm]["detected_at_1pct_false_alarm"]
                < lda["detector"][nm]["detected_at_1pct_false_alarm"])
    assert pct(rf["detector"]["radish"]["detected_at_1pct_false_alarm"]) == "93.9%"
    for nm in NAMES:   # the SVM is behind the LDA at 1 % false alarms on every species
        assert (svm["detector"][nm]["detected_at_1pct_false_alarm"]
                < lda["detector"][nm]["detected_at_1pct_false_alarm"])
    assert svm["n_errors"] == 276 and svm["millet_called_sorghum"] == 54
    assert pct(svm["detector"]["poplar"]["detected_at_1pct_false_alarm"]) == "14.6%"
    best = min(M[k]["n_errors"] for k in rows)
    assert best == 136
    assert S["separation"]["400-2400 nm"]["test_plants_unseen"]["n_errors"] == 23
    # "a fifth of its errors" in the qualified bullet
    assert 0.15 < 32 / lda["n_errors"] < 0.25
    assert "removes a fifth of its errors and none of the gap to 2400 nm" in README


def test_settings_recorded():
    s = R["models"]["settings"]
    assert "500 trees" in s["random forest"] and "C = 10" in s["rbf svm"]
    assert R["false_alarm_rates"] == [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]
    assert R["detection_levels"] == [0.90, 0.95, 0.99]
    assert re.search(r"Added 4 November 2026\.\*? Three things", README)
