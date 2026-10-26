# Star-catalog matching study

This repository records an exploratory research-practice exercise using **MATCH**, an existing C program by Michael W. Richmond for matching astronomical object catalogs under translation, rotation, scale, and distortion.

The work involved studying and testing this open-source tool during a research practice and considering its possible relevance to small-satellite tracking and space situational awareness. The repository is **not** a satellite-tracking system, and it does not contain orbit determination or a flight-ready implementation.

## Related later work at SPEL

A later, related development at the University of Chile's Space and Planetary Exploration Laboratory (SPEL) is **SOST (SPEL - Open Star Tracker)**, developed for the SUCHAI-2/3 and PlantSAT nanosatellite missions. The public [SPEL Star_Tracker repository](https://github.com/spel-uchile/Star_Tracker) explicitly lists **MATCH version 0.14** as one of the two open-source astronomy tools used as the basis of that star-tracker software (the other is Source Extractor). The doctoral thesis *Developing and testing an ultra-low-cost star tracker for attitude determinations of nanosatellites* (Samuel Tomás Gutiérrez Russell, 2023) documents this later work: [University of Chile repository record](https://repositorio.uchile.cl/handle/2250/194384).

This is a relevant technical connection, but the available public information does **not** establish that the SOST thesis or its implementation descended from this repository or from the specific experiments recorded here. Accordingly, this project is described as an earlier exploratory study of MATCH, not as a contribution to SOST. Any stronger claim would require confirmation from the thesis, its authors, or project records.

## Repository layout

- `analysis/coordinate_matching_exploration.ipynb` — a small notebook visualizing MATCH's bundled synthetic self-test catalogs and an example coordinate transformation.
- `third_party/match/` — the original MATCH source, documentation, build files, and test fixtures, retained as third-party software.
- `third_party/match/COPYING` — the upstream GNU General Public License, version 2.

## Run the notebook

Install Python with NumPy and Matplotlib, then open the notebook in Jupyter. Run it from either the repository root or the `analysis/` directory; it locates the bundled test data in either case.

The notebook is a visualization of **synthetic test fixtures bundled with MATCH**, not real satellite observations. Its example transformation coefficients are recorded values from an exploratory run; the notebook does not execute the MATCH program or estimate those coefficients itself.

## Attribution

MATCH was written by **Michael W. Richmond**. The upstream `AUTHORS`, `README`, documentation, source files, and license are retained under `third_party/match/`. The MATCH implementation is not original work by the repository owner. See the upstream license before redistributing or modifying that software.

## Scope and status

This is a cleaned record of an exploratory study, not a maintained standalone application. The notebook is intentionally limited to illustrating the catalog coordinates and an example transformation; it should not be interpreted as a validated satellite-tracking result.
