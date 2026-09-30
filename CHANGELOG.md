# Changelog

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
