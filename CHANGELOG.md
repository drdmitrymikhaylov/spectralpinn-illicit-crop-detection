# Changelog

## 2026-11-04

- The detector table's 1% false-alarm point read as a whole curve
  (`analysis/operating_points.py`, imports the September script; its linear
  rows reproduce `species_separation_vnir.json` to every digit). Millet goes
  53.1% → 82.3% → 97.9% at 1 / 2 / 5% false alarms and 13.5% at 0.1%; 90%
  detection costs 3.0% false alarms on millet and 0.6% or less on the four
  easy species; 99% costs 18.7% on pumpkin. A plant-level bootstrap (2 000
  resamples of the 128 plants) puts millet's 53.1% at [32.4%, 90.7%] and
  sorghum's 74.8% at [44.1%, 90.4%] — the two grass rows are one draw from
  16 plants each.
- The September threshold was placed on the same held-out leaves it was
  scored on. Placed on the training plants (inner out-of-fold posteriors) it
  delivers 0.5–1.3% false alarms and moves detection by at most 4.2 points
  (millet 57.3%); the in-sample threshold did not flatter the table.
- Two untuned non-linear models on the same plant-grouped folds, 400–1000 nm:
  a random forest makes 136 errors against the linear discriminant's 168
  (poplar → pumpkin 41 → 15; millet at 1% false alarms 53.1% → 78.1%, poplar
  93.7% → 86.9%) and leaves 42 of the 55 grass confusions; an RBF SVM makes
  276 and is behind the linear model on every species. The linear model in
  400–2400 nm leaves 23 — the gap between the windows is not a linear-model
  artefact. Two sentences of the September section qualified in place.
  `results/operating_points.json`, `tests/test_operating_points.py` (7
  tests, no dataset needed).

## 2026-10-19

- LICENSE now holds the standard CC BY 4.0 legal text, so GitHub recognises the licence (it showed "Other" before). The terms are unchanged; the note on what each licence covers moved to `NOTICE.md`.

## 2026-09-30

- First content beyond the description: a check on open leaf spectra of how far
  the 400–1000 nm window separates crop species (`analysis/species_separation_vnir.py`,
  numpy + scikit-learn; seven crops, 2 009 glasshouse leaves, EcoSIS data of
  Burnett et al. 2020, not redistributed). With test plants unseen the
  seven-way call is 91.6% right in 400–1000 nm against 98.9% in 400–2400 nm,
  168 errors against 23. At a 1% false-alarm rate detection runs from 99.4%
  (pepper) down to 53.1% (millet, whose AUC is still 0.987). In a linearly
  mixed pixel 90% detection takes a median target fraction of 0.40; at a tenth
  of the pixel no species pair reaches 90%, and three closely related pairs
  never reach it at all. A glasshouse leaf library calls field-grown pumpkin
  leaves correctly 58.5% of the time and above-canopy spectra 0 times in 62.
- New README section with a "what this does not show" list (none of the seven
  species is a target of the field work; linear mixing is the kindest model;
  no atmosphere, no flight), `results/species_separation_vnir.json`,
  `tests/test_readme_numbers.py` (5 tests, no dataset needed).
