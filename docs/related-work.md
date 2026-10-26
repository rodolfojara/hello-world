# Related work

## From exploratory study to later star-tracker research

This repository documents an exploratory study carried out during a research internship. The work investigated whether astronomical software for matching lists of stars could be useful in the context of small-satellite tracking.

The software studied was **MATCH**, developed by Michael W. Richmond. MATCH establishes relationships between lists of astronomical objects and can determine transformations between the coordinate systems of the lists.

The study explored a simple geometric approach: stellar positions could be projected mathematically and compared with positions obtained from a tracker, with the objective of investigating whether the correspondence could provide information about the spacecraft's position or pointing. The methods used at that stage were intentionally simple and exploratory.

### Later related work

The University of Chile doctoral thesis *Developing and testing an ultra-low-cost star tracker for attitude determinations of nanosatellites* later developed SOST, an open-source star tracker for nanosatellites. The thesis explicitly describes its Lost-In-Space algorithm as being based on the open-source astronomy software **Source Extractor and Match**.

The thesis explains that MATCH is used to establish the relationship between objects detected in an image and objects in a projected stellar-catalog segment. It also cites Michael W. Richmond as the author of MATCH.

This later work is relevant to the history of this repository, but there is no claim here that the author of this repository contributed code, ideas, or research to the doctoral thesis. The connection is simply that the same open-source astronomical tool studied in this earlier exploratory work was subsequently used as part of a much more developed star-tracker system.

## References

- Samuel Tomás Gutiérrez Russell, *Developing and testing an ultra-low-cost star tracker for attitude determinations of nanosatellites*, University of Chile, 2023.
- Michael W. Richmond, **MATCH — a program for matching star lists**.
- SOST / SPEL Open Star Tracker: https://github.com/spel-uchile/Star_Tracker
