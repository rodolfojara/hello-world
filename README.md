# Small-Satellite Tracking Study

Exploratory research work on astronomical star-list matching and geometric transformations, carried out during a research internship in the context of small-satellite tracking.

## About this repository

This repository preserves an early exploration of **MATCH**, an open-source astronomical program developed by Michael W. Richmond for establishing correspondences between lists of astronomical objects.

The original work investigated whether this type of software could be useful for associating observed stars with a reference catalog and for comparing observed tracker data with mathematically projected stellar positions. The approach was exploratory and used relatively simple geometric methods.

The purpose of this repository is to document that study and the technical direction it explored. It is **not** presented as an implementation of a complete star tracker.

## My contribution

My work consisted of studying and experimenting with existing open-source software, testing its behavior with sample star lists, exploring coordinate transformations, and considering possible improvements and applications to small-satellite tracking.

The MATCH software itself was not developed by me and remains third-party open-source software.

## Related later work

A later doctoral research project at the University of Chile developed SOST (SPEL Open Star Tracker), an open-source star tracker for nanosatellites. That work used established open-source astronomy software, including MATCH, as part of its Lost-In-Space attitude-determination algorithm.

This repository is **not part of the SOST development and does not claim any contribution to that doctoral thesis**. The later work is referenced only as related work: it demonstrates the relevance of the class of astronomical software and star-matching approach explored in this earlier study.

See [`docs/related-work.md`](docs/related-work.md) for more context.

## Repository status

This is a historical research record rather than a maintained software product. The repository is being cleaned up to distinguish third-party MATCH source code from the exploratory analysis carried out during the internship.
