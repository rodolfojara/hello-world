# Star-catalog matching study

This repository records an exploratory research-practice exercise using **MATCH**, an existing C program by Michael W. Richmond for matching astronomical object catalogs under translation, rotation, scale, and distortion.

The work involved studying and testing this open-source tool and considering its possible relevance to small-satellite tracking and space situational awareness. The repository is **not** a satellite-tracking system, and it does not contain orbit determination or a flight-ready implementation.

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
