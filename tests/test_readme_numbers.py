"""Pins the open-data section of README.md to results/species_separation_vnir.json.

Needs no dataset: it reads the page and the result file only.

    python -m pytest tests/test_readme_numbers.py
"""

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
README = (ROOT / "README.md").read_text(encoding="utf-8")
D = json.loads((ROOT / "results" / "species_separation_vnir.json").read_text(encoding="utf-8"))
NAMES = ["millet", "pepper", "poplar", "pumpkin", "radish", "sorghum", "sunflower"]


def on_page(s):
    assert s in README, f"not on the page: {s!r}"


def thousands(n):
    return f"{n:,}".replace(",", " ")


def test_the_data_are_what_the_page_says():
    d = D["data"]
    assert sorted(d["per_species"]) == NAMES
    assert sum(d["per_species"].values()) == d["n_glasshouse_leaf"]
    on_page(f"{thousands(d['n_glasshouse_leaf'])} glasshouse leaves from {d['plants']} "
            f"plants on {d['dates']} measurement dates")
    on_page(f"{d['band_step_nm']:.0f} nm bands")


def test_separation_table():
    for w, v in D["separation"].items():
        a, b = v["test_plants_unseen"], v["test_dates_unseen"]
        cm = a["confusion_rows_true"]
        assert a["n"] == sum(map(sum, cm)) == D["data"]["n_glasshouse_leaf"]
        assert a["n_errors"] == a["n"] - sum(cm[k][k] for k in range(7))
        on_page(f"| {w.replace('-', '–')} | {v['n_bands']} | {100 * a['accuracy']:.1f}% | "
                f"{a['n_errors']} | {100 * b['accuracy']:.1f}% |")
    vn = D["separation"]["400-1000 nm"]["test_plants_unseen"]
    full = D["separation"]["400-2400 nm"]["test_plants_unseen"]
    on_page(f"{vn['n_errors']} against {full['n_errors']}")
    assert 6.5 < vn["n_errors"] / full["n_errors"] < 7.5        # "by seven"
    cm = vn["confusion_rows_true"]
    mi, so, po, pu = (NAMES.index(n) for n in ("millet", "sorghum", "poplar", "pumpkin"))
    on_page(f"Of the {vn['n_errors']}, {cm[mi][so] + cm[so][mi]} are the two C4 grasses")
    on_page(f"millet called sorghum {cm[mi][so]} times in {sum(cm[mi])}, sorghum called "
            f"millet {cm[so][mi]} times in {sum(cm[so])}")
    on_page(f"{cm[po][pu]} are poplar called pumpkin")


def test_detector_table():
    det = D["detector"]["400-1000 nm"]
    for s, x in det.items():
        assert x["missed_of_target"] == round(
            (1 - x["detected_at_1pct_false_alarm"]) * x["n_target"])
        on_page(f"| {s} | {x['auc']:.3f} | {100 * x['detected_at_1pct_false_alarm']:.1f}% | "
                f"{x['missed_of_target']} of {x['n_target']} |")
    worst_auc = min(x["auc"] for x in det.values())
    on_page(f"Every AUC is {worst_auc:.3f} or better")
    m = det["millet"]
    assert m["auc"] == worst_auc
    on_page(f"misses {m['missed_of_target']} leaves of {m['n_target']}")
    full = min(x["detected_at_1pct_false_alarm"] for x in D["detector"]["400-2400 nm"].values())
    on_page(f"every species is detected at {100 * full:.1f}% or better")


def test_mixed_pixel_table():
    m = D["mixed_pixel"]
    S, fr = m["summary"], m["fractions"]
    assert m["n_pairs"] == 42 and m["false_alarm_rate"] == 0.01
    assert fr == [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
    for label, key in (("detected, median pair", "median"), ("detected, worst pair", "min"),
                       ("detected, best pair", "max")):
        on_page(f"| {label} | " + " | ".join(f"{100 * S[str(f)][key]:.0f}%" for f in fr) + " |")
    on_page("| pairs at 90% or better, of 42 | "
            + " | ".join(str(S[str(f)]["pairs_at_or_above_90pct"]) for f in fr) + " |")
    assert S["0.1"]["pairs_at_or_above_90pct"] == 0
    f90 = S["fraction_for_90pct"]
    reach = {k: v["fraction_for_90pct"] for k, v in m["pairs"].items()
             if v["fraction_for_90pct"] is not None}
    never = {k: v["detected"]["1.0"] for k, v in m["pairs"].items()
             if v["fraction_for_90pct"] is None}
    assert len(reach) == 39 and len(never) == f90["pairs_never_reaching_it"] == 3
    lo, hi = min(reach, key=reach.get), max(reach, key=reach.get)
    on_page(f"a median fraction of {f90['median']:.2f} over the {len(reach)} pairs that get "
            f"there, from {reach[lo]:.2f} ({lo.replace(' in ', ' among ')}) to "
            f"{reach[hi]:.2f} ({hi.replace(' in ', ' among ')})")
    for k in ("millet in sorghum", "sorghum in millet", "poplar in pumpkin"):
        assert k in never
        on_page(f"{k.replace(' in ', ' among ')} ({100 * never[k]:.0f}%)")


def test_new_venue_table():
    v = D["new_venue"]
    cm = D["separation"]["400-1000 nm"]["test_plants_unseen"]["confusion_rows_true"]
    pu = NAMES.index("pumpkin")
    assert abs(cm[pu][pu] / sum(cm[pu]) - v["glasshouse_recall_test_plants_unseen"]) < 1e-12
    on_page(f"| glasshouse leaves, plants unseen | "
            f"{100 * v['glasshouse_recall_test_plants_unseen']:.1f}% | "
            f"{cm[pu][pu]} of {sum(cm[pu])} |")
    leaf, can = v["field_leaf"], v["field_canopy"]
    on_page(f"| field leaves, same leaf clip | {100 * leaf['rate']:.1f}% | "
            f"{leaf['called_target']} of {leaf['n']} |")
    on_page(f"| field, from above the canopy | {100 * can['rate']:.1f}% | "
            f"{can['called_target']} of {can['n']} |")
    assert 0.38 < 1 - leaf["rate"] < 0.45                      # "four leaves in ten"
    c = can["called"]
    assert sum(c.values()) == can["n"] and c["pumpkin"] == 0
    on_page(f"called pepper ({c['pepper']}), radish ({c['radish']}), "
            f"sunflower ({c['sunflower']}) and sorghum ({c['sorghum']})")
    assert c["millet"] == c["poplar"] == 0
