# GloriPaka packaging roadmap

Turning the Strategy-B research workspace into a clean, reproducible, standalone Python
package. Checked items are complete; the rest are planned. Built minimal-first, then
expanded.

## Phase 1 - Minimal viable package (complete)
- [x] Clean repo tree; package renamed from `strategyb` to `gloripaka`.
- [x] `config.py` decoupled from the parent repository (self-contained, repo-relative paths).
- [x] `pyproject.toml` for uv (metadata, dependencies, optional extras, ruff + mypy + pytest, Google docstring convention).
- [x] `LICENSE` (MIT), `CITATION.cff`, `NOTICE`, `.gitignore`, `.gitattributes`, `py.typed`.
- [x] `README.md` with logo, badges, installation, usage, module table, and citation.
- [x] `uv sync` builds and installs `gloripaka`; imports and a functional smoke test pass.
- [ ] Push to `GloriPaka` on GitHub.

## Phase 2 - End-to-end pipeline (sequence to prediction)
- [ ] Port featurization into `src/gloripaka/featurize/`; deduplicate against the World-1 copy; keep the DeepSP weights and license.
- [ ] `src/gloripaka/structure.py`: NanoBodyBuilder2 structure generation from sequence.
- [ ] `src/gloripaka/pipeline.py`: high-level `GloriPakaPredictor` (sequences in, HIC and confidence out).
- [ ] Command-line entry point (`gloripaka predict input.csv`).
- [ ] Curate `experiments/` (keep the reproducible experiments; each runs against in-repo data).
- [ ] `scripts/reproduce.py` to regenerate the headline result and key figures.

## Phase 3 - Data and Zenodo
- [ ] Consolidate the small canonical CSVs into `data/` (inputs, labels, and the feature tables actually used).
- [ ] Keep a small in-repo sample (the 34 known VHH) under `data/examples/` for tests and examples.
- [ ] Stage the heavy artifacts (PDBs, full feature tables) in the separate `zenodo_gloripaka/` folder; write deposition metadata; upload; add the DOI and a `download_data.py` helper.
- [ ] `data/README.md` documenting provenance and how to fetch.

## Phase 4 - Literature and dataset construction
- [ ] Curate the literature notes into `docs/literature.md` with citations.
- [ ] `docs/datasets.md`: how each dataset (the World-1 pool, native VHH, scale bridges, and the ARRAS test set) was constructed, with DOIs.
- [ ] Decide manuscript inclusion carefully in view of figure and publisher licensing; likely link out rather than commit.

## Phase 5 - Tests, typing, documentation, CI
- [ ] `tests/`: smoke and unit tests (metrics, physics composite, within-group quantiles, clone-family clustering, model fit and predict on a toy set).
- [ ] Docstring and typing pass; `ruff check` and `mypy src/` clean.
- [ ] Documentation site (mkdocs-material) published to GitHub Pages.
- [ ] GitHub Actions running ruff, mypy, and pytest on push and pull request.
- [ ] `CONTRIBUTING.md` and issue templates.

## Phase 6 - Release
- [ ] Tag `v0.1.0` and create a GitHub release; finalize the Zenodo, docs, and CI badges; add README figures and cross-links.

---
The source of truth for the science is the research lab journal
(`notes/lab_journal.md` in the source workspace). This roadmap tracks only the
packaging and reproducibility effort.
