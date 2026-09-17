# The multilinear gauge-invariant operator basis for five massless photons

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22819347.svg)](https://doi.org/10.5281/zenodo.22819347)

**Status: preprint, not peer-reviewed.** See the disclaimer on the first page of the paper.
Feedback, corrections, and criticism are genuinely welcome — please open an issue.

## What this is

Extends the S-matrix classification program of Henning, Lu, Melia & Murayama
(arXiv:1706.08520) and Chowdhury, Gadde, Gopalka, Halder, Janagal & Minwalla (arXiv:1910.14392)
to the spinning (photon) case at `n=5` particles:

- The `S_5` representation on the space of independent Mandelstam invariants (irreducible,
  dimension 5) and its complete invariant ring (cross-validated exactly against the published
  scalar-case result).
- The gauge-invariant photon polarization data structure at `n=5` (two physical parameters per
  particle, not one as at `n=4`), with closed-form expressions.
- Four families of gauge-invariant building blocks, including a new one (the "triple
  field-strength sandwich").
- A proved parity theorem, and a classification of multilinear five-photon invariants through
  mass dimension 11.

This repository does **not** address the Classical Regge Growth conjecture at `n≥5` — that is a
separate, harder open question, treated in a companion note (in preparation).

## Structure

```
paper/    — the LaTeX source and compiled PDF (start here)
scripts/  — the Python/Singular scripts that produced every number cited in the paper
notes/    — session notes documenting how each result was derived, including bugs found and
            fixed along the way (referenced from the paper's verification-methodology section)
```

## Reproducing the results

Every script in `scripts/` is self-contained and runnable with a standard Python 3 + `sympy` +
`mpmath` install; a few require Singular (via WSL on Windows) for the invariant-ring computation
in `s5_invariant_ring.sing`. There is no single "run everything" entry point — each script
corresponds to a specific claim in the paper; see the comments at the top of each file for what it
verifies.

## License

Code (`scripts/`) is MIT-licensed. Text and the manuscript (`paper/`, `notes/`) are licensed
under CC-BY 4.0. See `LICENSE`.

## Provenance

This repository was extracted from a larger, private working repository where this result was
developed alongside other exploratory work; the history here starts fresh from that point rather
than replaying every intermediate commit.
