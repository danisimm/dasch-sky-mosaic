# Archived: SExtractor + SCAMP single-image WCS solver

Shelved in favour of the astroalign two-image route (JPG + FITS mosaic, see
`pipeline.py` / `wcs_solve.py` on `main`). Work happened Jul 10–14 2026.

## Idea
Solve the plate JPG's WCS from the JPG alone (one download per plate instead
of two): rough initial WCS from the exposures catalog -> SExtractor detections
-> SCAMP fit against a local Gaia DR3 G<14 catalog, trying 4 rotation
hypotheses (0/90/180/270) and keeping the lowest reduced chi2 with >= 20
matched detections.

## Files on this branch
- `src/dasch_sky_mosaic/wcs_solve_scamp.py` – solver + CLI (runs in WSL only;
  needs `source-extractor` and `scamp` installed)
- `src/dasch_sky_mosaic/scamp_prep.py` – initial-guess WCS (empirical scale
  correction 1.19, always mirror-flipped) and JPG->FITS conversion
- `src/dasch_sky_mosaic/scamp_config/` – SExtractor/SCAMP config
- `dev/gaia_catalog.py` – downloads full-sky Gaia DR3 G<14 in dec bands (astroquery)
- `dev/gaia_to_ldac.py` – converts it to SCAMP's LDAC reference format
- `dev/build_single_exposure_catalog.py` – `exposures.csv` -> `exposures_single.csv`
  (the initial-guess input)
- `dev/gather_stats.py` – compared the initial guess against astroalign solutions
  (-> `scamp_seed_analysis.csv`)
- `dev/run_scamp_batch.py` – run the solver over every cached JPG
- `compare_sources.py` – sep vs SExtractor detection comparison (debugging)
- `pyproject.toml` – adds `astroquery`

## Where it got to
Ran over 23 cached plates; 13 produced a solution. Plate centres vs astroalign:

| plate | n_matched | reduced chi2 | centre diff vs astroalign |
|---|---|---|---|
| b02312 | 550 | 0.14 | 0.53' |
| b31258 | 136 | 0.81 | 0.43' |
| dnb06698 | 32 | 0.64 | 2.91' |
| dnb06700 | 59 | 0.11 | 0.97' |
| dnb06706 | 26 | 0.07 | 0.48' |
| dnb06709 | 77 | 0.10 | (no astroalign) |
| dsb02842 | 90 | 0.07 | 1.07' |
| dsb02843 | 30 | 0.42 | 1.01' |
| dsb02887 | 25 | 0.42 | 0.89' |
| h01481 | 21 | 114 (bad) | (no astroalign) |
| j01752 | 35 | 5.1, all 4 rotations "solved" | (no astroalign) |
| kd00186 | 27 | 0.03 | 5.76' |
| rh14460 | 143 | 0.07 | 0.54' |

Open problem when shelved: SCAMP was slow. An strace of one run showed ~99%
of wall time in `futex` (thread waits), not I/O.

## Data (not in git)
Moved to `data/archive/scamp/` in the working copy (gitignored):
- `reference_catalogs/` – Gaia DR3 G<14 raw bands + merged + LDAC (~2.2 GB;
  slow to re-download)
- `scamp_solve_cache/` – per-plate intermediates (~12 GB; regenerable, safe to delete)
- `wcs_scamp/` – solver output JSON
- `scamp_seed_analysis.csv`, `exposures_single.csv`
- `debug/` – `scamp.xml`, `scamp_rot0.perf`, `strace_rot0.txt`

To revive: `git switch archive/scamp-single-image`, then move the data back to
`data/reference_catalogs/` and `dev/exposures_single.csv`.
