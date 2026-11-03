# Perfect play with spatial quantum strategies is Σ⁰₂-complete

Seth Douglas and Nidhal Mghirbi. Version 1.1.0, November 2026. Concept DOI (all versions): [10.5281/zenodo.23910986](https://doi.org/10.5281/zenodo.23910986). Version 1.0.0: [10.5281/zenodo.23910987](https://doi.org/10.5281/zenodo.23910987). Version 1.1.0: [10.5281/zenodo.23929983](https://doi.org/10.5281/zenodo.23929983).

**This is a preprint. It has not been peer reviewed.**

## Abstract

When does a finite nonlocal game admit perfect spatial quantum play? We prove that this decision problem is complete for the second existential level of the arithmetical hierarchy. The same classification holds when a perfect strategy must have finite marginal entropy. Hardness persists among games whose finite-dimensional value is one and which have no finite-dimensional perfect strategy. A finite clock-phase order supplies the existential witness: normality forces a finite order on a nonzero sector, while the strategy chooses the order.

For the upper bound, we normalize the largest Schmidt coefficient and use compact operator sandwiches. Fixing a lower bound on the largest Schmidt probability gives a co-recursively enumerable stratum. No computable bound selects a sufficient stratum for every game. Exact attainment of a rational Bell bound, of an unspecified supremum, and of a specified scalar value is also Sigma^0_2-complete. The hardness proof extends the modular-shift and rational-row toolkit of our earlier entropy-game paper to four calibrated bridges and a deterministic computation on rooted finite squares.

## Contents and reproduction

- `paper/main.pdf`: manuscript.
- `paper/`: LaTeX source, bibliography, plain-text abstract, source manifest, source archive and build driver. Run `python -B build.py` there with pdflatex and BibTeX installed.
- `verification/`: exact and numerical finite checks, recorded outputs and reproduction commands. Start with `verification/README.txt`. `verification/SHA256SUMS` records the bundled verification files.
- The only compiled game is the minimal compiler instance: four bridges, one modular isometry and no transported labels, with 748 rational rows and 27,991 questions per party. The full reduction game was not compiled.

Computational checks support the written proofs. Finite evaluations alone do not prove the general compiler, the infinite-dimensional compactness argument or the hardness reduction.

## Changes in version 1.1.0

- No theorem statement has changed.
- Section 4 (the compiler) adds a one-page overview, "From a presentation to a game", and writes out in full the infinite-dimensional steps in the proofs of Lemmas 4.1 to 4.4: normal states and purifications, trace-class densities, supports, the dense-range cancellations, and the exact operator identities used for soundness and completeness.
- Section 10.3, now "The von Neumann type ladder", replaces the questions posed in version 1.0.0 by the results of three papers by the authors: a finite game that type III algebras win perfectly and semifinite algebras cannot, a game that normal hyperfinite strategies approach but never win, and the arithmetical classification of perfect play by algebra type.
- Section 10.4 adds that perfect spatial play attainable only with infinite entropy is D₂(Σ⁰₂)-complete (Corollary 9.2 of the type-ladder paper).
- The paper, this README and `verification/README.txt` state that the compiled 27,991-question game is the minimal instance and that the full reduction game was not compiled.

## Related work

- *Perfect play in a finite nonlocal game requires infinite entanglement entropy*, whose toolkit this paper extends: https://doi.org/10.5281/zenodo.23902656
- *Finite games that read the type of a von Neumann algebra*: https://doi.org/10.5281/zenodo.23787124
- *A finite game that normal hyperfinite strategies approach but never win*: https://doi.org/10.5281/zenodo.23816141
- *Perfect play has different arithmetical complexity in semifinite and type III algebras*: https://doi.org/10.5281/zenodo.23925092

## Citation and rights

See `CITATION.cff` for citation metadata.
Paper and documentation: CC BY 4.0. Scripts: MIT or Apache 2.0, at your option. See `RIGHTS.md`, the license files and `NOTICE`.
