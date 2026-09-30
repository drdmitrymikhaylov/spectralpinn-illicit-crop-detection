# Spectral PINN — Detecting Illicit Crops from UAV Spectra

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

Every AUC is 0.987 or better. At the operating point a survey needs, the millet detector with that AUC misses 45 leaves of 96. An AUC near 1 is not a detection rate. With the range extended to 2400 nm every species is detected at 98.7% or better.

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
- The detector is a linear discriminant. It is a baseline, not the radiative-transfer model this page describes.

## Where this stands

This is deployed work with a United Nations agency, and a patent has been filed. This page is the public description — the operational pipeline and the survey data are not public. The source is available under NDA where appropriate. What the repository does hold is the open-data reference calculation above (`analysis/`), its result file (`results/`) and the tests that pin this page to it (`tests/`); see `CHANGELOG.md`.

## Licence

Documentation, figures, result files and the reference calculation in `analysis/` in this repository: CC BY 4.0. Source code of the operational pipeline is held in a private repository, all rights reserved, and is available under NDA. The open spectra belong to their authors and are not redistributed here: A. C. Burnett, S. P. Serbin, K. J. Davidson, K. S. Ely, A. Rogers. 2020. *Hyperspectral leaf reflectance, biochemistry, and physiology of droughted and watered crops.* Data set, EcoSIS.
