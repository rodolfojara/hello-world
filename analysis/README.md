# Exploratory analysis

This directory contains the author's exploratory analysis from the small-satellite tracking research practice.

## Notebook

- [`coordinate-matching-exploration.ipynb`](coordinate-matching-exploration.ipynb): reads and plots two sample star lists, then applies a simple affine coordinate transformation to one list for visual comparison.

The sample files are part of the third-party MATCH distribution and are located in `../third_party/MATCH/testdata/`.

## Dependencies and limitations

The notebook uses Python, Matplotlib, and Astropy (`astropy.io.ascii`). It is a historical exploratory notebook, not a complete star-tracker implementation. The transformation coefficients are recorded in the notebook; the notebook does not document a full estimation pipeline from tracker measurements to spacecraft attitude or position.
