# Spectral PINN — Detecting Illicit Crops from UAV Spectra

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23529573.svg)](https://doi.org/10.5281/zenodo.23529573)

## The field work behind it

This is field technology. UAV spectrometers flown over hidden cannabis and ephedra plantations in rough terrain with a United Nations agency. Published in 2025 (*UAV mounted spectrometry camera solution for detecting illicit cannabis and ephedra plantations*, Springer), patent application filed, deployed. This page explains how a spectrometer in the air picks out an illicit crop from everything growing around it.

## Problem

Narcotic crops grow inside legal fields, in places nobody can survey on foot. A 400–1000 nm spectrometer flown at 120 m separates them from the crops around them. The reflectance signature is a property of the leaf, not of the field boundary.

## Approach

Radiative transfer through the canopy. What a spectrometer measures is scattering and absorption by pigments, water and cell structure, integrated over the leaf stack and the soil behind it. Species separation then turns on two questions. Which absorption features survive that integration at flight altitude, and how far can an atmospheric correction be trusted between calibration flights?

### Related work in this line

- UAV mounted spectrometry camera solution for detecting illicit plantations — *Discover Applied Sciences*.
- Adaptive Calibration Cycle — keeping a hyperspectral system honest between laboratory calibrations.

## How far does 400–1000 nm separate species? A check on open data

*Added 30 September 2026.* The survey data behind this page is not public, so a reader could check nothing on it. This section puts numbers on the first question under "Approach" using open leaf spectra of seven ordinary crops: pumpkin, pepper, poplar, radish, sunflower, sorghum and foxtail millet (Burnett, Serbin, Davidson, Ely and Rogers, 2020, EcoSIS). That is 2 009 glasshouse leaves from 128 plants on 65 measurement dates, one spectroradiometer, 4 nm bands. None of the seven is a target of the field work. They stand in for "one species among others". `analysis/species_separation_vnir.py` computes everything below, `results/species_separation_vnir.json` holds it, and `tests/test_readme_numbers.py` pins this section to the file.

**Seven species, test plants never seen in training.** A linear discriminant on brightness-normalised spectra:

| window | bands | accuracy | errors of 2 009 | accuracy, whole dates held out |
|---|---|---|---|---|
| 400–1000 nm | 151 | 91.6% | 168 | 91.1% |
| 400–2400 nm | 501 | 98.9% | 23 | 98.6% |

Stopping at 1000 nm multiplies the errors by seven, 168 against 23. Holding out whole measurement dates changes little, so the model is not reading the day. The errors are not spread evenly. Of the 168, 55 are the two C4 grasses taken for each other (millet called sorghum 23 times in 96, sorghum called millet 32 times in 151) and 41 are poplar called pumpkin.

**As a detector at a 1% false-alarm rate.** Each species in turn against the other six, 400–1000 nm:

| target | AUC | detected at 1% false alarms | missed |
|---|---|---|---|
| pepper | 1.000 | 99.4% | 3 of 498 |
| radish | 0.999 | 96.7% | 7 of 214 |
| sunflower | 0.998 | 96.5% | 6 of 172 |
| poplar | 0.998 | 93.7% | 26 of 411 |
| pumpkin | 0.988 | 86.7% | 62 of 467 |
| sorghum | 0.991 | 74.8% | 38 of 151 |
| millet | 0.987 | 53.1% | 45 of 96 |

Every AUC is 0.987 or better. At the operating point a survey needs, the millet detector with that AUC misses 45 leaves of 96. An AUC near 1 is not a detection rate. With the range extended to 2400 nm every species is detected at 98.7% or better. *(Added 4 November 2026: the 1% column is one point on a steep curve — at 2% false alarms millet is at 82.3%, at 5% at 97.9% — and, drawn from 16 plants, it carries a plant-level interval of 32.4% to 90.7%. The whole curve is in "Operating points" below.)*

**A mixed pixel.** The field problem is a target crop planted among a host crop, so a pixel is part target. Here target and host leaf reflectance are mixed linearly by area fraction, the detector is trained on pure leaves of other plants, the threshold gives 1% false alarms on pure host, and all 42 ordered pairs of species are run, 400–1000 nm:

| target fraction of the pixel | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.7 | 1.0 |
|---|---|---|---|---|---|---|---|
| detected, median pair | 14% | 37% | 71% | 86% | 97% | 100% | 100% |
| detected, worst pair | 2% | 3% | 4% | 5% | 9% | 29% | 71% |
| detected, best pair | 52% | 94% | 100% | 100% | 100% | 100% | 100% |
| pairs at 90% or better, of 42 | 0 | 2 | 10 | 20 | 29 | 33 | 39 |

With a tenth of the pixel no pair reaches 90%. Ninety per cent detection takes a median fraction of 0.40 over the 39 pairs that get there, from 0.19 (sorghum among sunflower) to 0.87 (pumpkin among poplar). Three pairs never get there, even with a pure pixel: millet among sorghum (71%), sorghum among millet (86%) and poplar among pumpkin (89%). A close relative as the host is a different problem from an unrelated crop as the host.

**The same species somewhere else.** The glasshouse model, asked about pumpkin grown in a field:

| pumpkin spectra | called pumpkin | |
|---|---|---|
| glasshouse leaves, plants unseen | 94.6% | 442 of 467 |
| field leaves, same leaf clip | 58.5% | 196 of 335 |
| field, from above the canopy | 0.0% | 0 of 62 |

A library built on glasshouse leaves loses four pumpkin leaves in ten once the plant grows outdoors, and recognises none of the 62 spectra taken from above. Those are called pepper (31), radish (17), sunflower (13) and sorghum (1). The sentence under "Problem", that the signature is a property of the leaf, needs this qualification: a spectrometer in the air sees a canopy, and on these data a leaf library did not carry to it. The canopy needs its own calibration.

**What this does not show.**

- None of the seven species is cannabis or ephedra, and no survey spectrum is used. How separable the real targets are from their real hosts is not answered here.
- Linear mixing of leaf reflectance is the kindest mixing model there is. A real pixel adds soil, shadow, leaf angle and scattering between leaves, and all of them move the numbers the wrong way. The fractions in the table are the least a field pixel needs.
- One instrument on a leaf clip with its own lamp. No atmosphere, no flight. The second question under "Approach", how far an atmospheric correction can be trusted between calibration flights, is not touched.
- The species were not randomised across growing conditions. Each came from its own pots and its own schedule. Holding out dates guards against reading the day, not against reading the experiment, so the separation table is an upper bound for that reason too.
- The detector is a linear discriminant. It is a baseline, not the radiative-transfer model this page describes. *(Added 4 November 2026: a random forest on the same folds removes a fifth of its errors and none of the gap to 2400 nm — see below.)*

## Operating points, a threshold set without the test leaves, and two non-linear models

*Added 4 November 2026.* Three things the section above left open, on the same 2 009 leaves, the same quality gate and the same plant-grouped folds. `analysis/operating_points.py` computes them (it imports the script above, and its linear-discriminant rows reproduce that script's numbers to every digit), `results/operating_points.json` holds them, `tests/test_operating_points.py` pins this section to the file.

**The 1% column is one point on a steep curve.** Detection of each species against the other six, 400–1000 nm, at six false-alarm rates; the last column resamples the 128 plants 2 000 times and reads the 1% column again each time:

| target | 0.1% | 0.5% | 1% | 2% | 5% | 10% | 1% column, plants resampled |
|---|---|---|---|---|---|---|---|
| pepper | 98.8% | 99.4% | 99.4% | 99.6% | 99.8% | 99.8% | [98.7%, 100.0%] |
| radish | 93.5% | 95.8% | 96.7% | 98.6% | 99.5% | 100.0% | [93.0%, 99.6%] |
| sunflower | 62.2% | 94.2% | 96.5% | 98.8% | 100.0% | 100.0% | [92.8%, 99.5%] |
| poplar | 77.9% | 88.6% | 93.7% | 98.1% | 99.8% | 100.0% | [85.4%, 98.9%] |
| pumpkin | 69.4% | 81.2% | 86.7% | 90.8% | 94.4% | 98.1% | [80.1%, 93.3%] |
| sorghum | 43.7% | 64.9% | 74.8% | 82.1% | 97.4% | 100.0% | [44.1%, 90.4%] |
| millet | 13.5% | 39.6% | 53.1% | 82.3% | 97.9% | 100.0% | [32.4%, 90.7%] |

Doubling the false-alarm budget from 1% to 2% takes millet from 53.1% to 82.3%; at 5% every species is above 94%, at 10% above 98%. Going the other way is expensive: a survey that can afford one false alarm in a thousand pixels keeps 98.8% of pepper and 13.5% of millet, and sunflower falls from 94.2% to 62.2% between 0.5% and 0.1%. The interval column is the part the 30 September table did not say: 53.1% for millet is one draw from 16 plants, and resampling plants puts it anywhere from a third to nine in ten. Sorghum's 74.8% reads [44.1%, 90.4%] the same way. The intervals of the four easy species stay above 85%; pumpkin's, from 24 plants, is [80.1%, 93.3%].

**What a survey has to pay for a detection rate, and whether the threshold was honest.** The first three columns read the same curves the other way: the false-alarm rate at which 90%, 95% and 99% of the target leaves are caught. The last two columns repeat the 1% operating point with the threshold set on the training plants only (inner out-of-fold posteriors), so the held-out leaves are scored and never used to place it:

| target | false alarms for 90% | for 95% | for 99% | threshold from training plants: false alarms delivered | detected |
|---|---|---|---|---|---|
| pepper | 0.0% | 0.0% | 0.3% | 0.5% (7 of 1 511) | 99.4% |
| radish | 0.1% | 0.3% | 2.6% | 0.9% (16 of 1 795) | 96.3% |
| sunflower | 0.3% | 0.7% | 1.9% | 1.0% (18 of 1 837) | 96.5% |
| poplar | 0.6% | 1.2% | 2.5% | 0.9% (15 of 1 598) | 92.0% |
| pumpkin | 1.8% | 5.3% | 18.7% | 1.3% (20 of 1 542) | 86.7% |
| sorghum | 2.9% | 4.1% | 5.3% | 0.9% (17 of 1 858) | 73.5% |
| millet | 3.0% | 4.8% | 5.2% | 1.0% (19 of 1 913) | 57.3% |

Ninety per cent detection costs 3.0% false alarms on millet, 2.9% on sorghum, 1.8% on pumpkin and 0.6% or less on the other four; 99% costs 18.7% on pumpkin, 5.2–5.3% on the two grasses and 2.6% or less on the rest. The threshold in the 30 September table was placed on the same held-out leaves it was scored on. Placed on the training plants instead, it delivers 0.5–1.3% false alarms (pumpkin the highest, 20 leaves of 1 542) and moves detection by at most 4.2 points — millet 53.1% → 57.3%, poplar 93.7% → 92.0%, sorghum 74.8% → 73.5%, radish 96.7% → 96.3%, the other three unchanged. The in-sample threshold did not flatter the table. In 400–2400 nm the same honest threshold delivers 0.6–1.1% false alarms and 98.9–100% detection (millet 95 of 96).

**Two non-linear models on the same folds.** The last bullet above called the linear discriminant a baseline. A random forest (500 trees) and an RBF support-vector machine (C = 10), both untuned, on the same SNV spectra and the same five plant-grouped folds, 400–1000 nm:

| model | accuracy | errors of 2 009 | millet ↔ sorghum | poplar → pumpkin | millet at 1% | sorghum at 1% | pumpkin at 1% | poplar at 1% |
|---|---|---|---|---|---|---|---|---|
| linear discriminant | 91.6% | 168 | 55 | 41 | 53.1% | 74.8% | 86.7% | 93.7% |
| random forest | 93.2% | 136 | 42 | 15 | 78.1% | 76.2% | 82.2% | 86.9% |
| RBF SVM | 86.3% | 276 | 71 | 39 | 47.9% | 68.9% | 61.0% | 14.6% |

The forest removes 32 of the 168 errors, 26 of them poplar-called-pumpkin (41 → 15), and takes the millet detector at 1% false alarms from 53.1% to 78.1% — but it gives back detection on poplar (93.7% → 86.9%), pumpkin (86.7% → 82.2%) and radish (96.7% → 93.9%), and it leaves 42 of the 55 grass confusions in place. The support-vector machine at one untuned setting is behind the linear model on every species at 1% false alarms (poplar 14.6%) and makes 276 errors, 54 of them millet called sorghum. Against the full range the comparison that matters: the best model in 400–1000 nm leaves 136 errors where the linear one in 400–2400 nm leaves 23. The gap between the two windows is not a linear-model artefact; the information is not in the bands.

**What this does not show.**

- One setting per non-linear model. A tuned SVM would very likely do better than 86.3%, and a forest with tuned depth and features might do better than 93.2%; neither tuning was run, and either would need a nested split so that the tuning does not see the test plants.
- The plant resampling resamples the scoring, not the training: the folds and the models are those of the 30 September run, so the intervals say how much the 1% column depends on which plants are held out, not how much a retrained model would move.
- The honest threshold is honest about the glasshouse. The 1.0% it delivers on millet is on glasshouse leaves of unseen plants; what the same threshold delivers on field leaves is the subject of the pumpkin table above (58.5%), not of this one.
- False alarms here are leaves of the other six species. A survey's non-targets are soil, shadow, water and every plant not in the library, and the false-alarm rate against those is not measured by anything on this page.

## Where this stands

This is deployed work with a United Nations agency, and a patent has been filed. This page is the public description — the operational pipeline and the survey data are not public. The source is available under NDA where appropriate. What the repository does hold is the open-data reference calculation above (`analysis/`), its result file (`results/`) and the tests that pin this page to it (`tests/`); see `CHANGELOG.md`.

## Licence

Documentation, figures, result files and the reference calculation in `analysis/` in this repository: CC BY 4.0. Source code of the operational pipeline is held in a private repository, all rights reserved, and is available under NDA. The open spectra belong to their authors and are not redistributed here: A. C. Burnett, S. P. Serbin, K. J. Davidson, K. S. Ely, A. Rogers. 2020. *Hyperspectral leaf reflectance, biochemistry, and physiology of droughted and watered crops.* Data set, EcoSIS.
