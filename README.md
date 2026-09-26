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

## Where this stands

This is deployed work with a United Nations agency, and a patent has been filed. This page is the public description — the operational pipeline and the survey data are not public. The source is available under NDA where appropriate.

## Licence

Documentation, figures and result files in this repository: CC BY 4.0. Source code is held in a private repository, all rights reserved, and is available under NDA.
